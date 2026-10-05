// Full-history p=0 Kalman filter/backward sampler with MATLAB-supplied noise.
// B = dsc_mean_mex(Y,h,P,mask,V,prior_mean,prior_covariance,noise)
#include "mex.h"
#include "lapack.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
struct Failure : std::runtime_error {
    const char* id;
    Failure(const char* identifier,const char* message):std::runtime_error(message),id(identifier) {}
};
void require(bool condition,const char* message) {
    if (!condition) throw Failure("dsc:MeanInput",message);
}
const double* values(const mxArray* a) {
    require(mxIsDouble(a)&&!mxIsComplex(a)&&!mxIsSparse(a),"Numeric inputs must be real full doubles.");
    return mxGetDoubles(a);
}
void symmetric(std::vector<double>& S,mwSize m) {
    for (mwSize j=0;j<m;++j) for (mwSize i=j+1;i<m;++i) {
        const double v=.5*(S[i+m*j]+S[j+m*i]); S[i+m*j]=v; S[j+m*i]=v;
    }
}
std::vector<double> right_solve(const std::vector<double>& A,
    std::vector<double> F,mwSize m,mwSize n) {
    std::vector<double> rhs(n*m);
    for (mwSize j=0;j<m;++j) for (mwSize i=0;i<n;++i) rhs[i+n*j]=A[j+m*i];
    const char lower='L'; const ptrdiff_t dimension=n, columns=m; ptrdiff_t info=0;
    dposv(&lower,&dimension,&columns,F.data(),&dimension,rhs.data(),&dimension,&info);
    if (info!=0) throw Failure("dsc:StateCovariance","Kalman prediction covariance is not positive definite.");
    std::vector<double> gain(m*n);
    for (mwSize j=0;j<n;++j) for (mwSize i=0;i<m;++i) gain[i+m*j]=rhs[j+n*i];
    return gain;
}
void normal_draw(std::vector<double> S,const std::vector<double>& mu,
    const double* noise,double* output,mwSize t,mwSize T,mwSize m) {
    symmetric(S,m); const char lower='L',vectors='V'; const ptrdiff_t n=m; ptrdiff_t info=0;
    std::vector<double> original=S;
    dpotrf(&lower,&n,S.data(),&n,&info);
    if (info==0) {
        for (mwSize i=0;i<m;++i) {
            double value=mu[i];
            for (mwSize k=0;k<=i;++k) value+=S[i+m*k]*noise[k];
            output[t+T*i]=value;
        }
        return;
    }
    S=original; std::vector<double> eigenvalues(m),work(std::max<mwSize>(1,3*m));
    const ptrdiff_t workspace=work.size();
    dsyev(&vectors,&lower,&n,S.data(),&n,eigenvalues.data(),work.data(),&workspace,&info);
    if (info!=0) throw Failure("dsc:StateCovariance","State covariance eigensolver failed.");
    double scale=1;
    for (double eigenvalue:eigenvalues) scale=std::max(scale,std::abs(eigenvalue));
    for (double eigenvalue:eigenvalues)
        if (eigenvalue < -1e-10*scale)
            throw Failure("dsc:StateCovariance","State covariance has a material negative eigenvalue.");
    for (mwSize i=0;i<m;++i) {
        double value=mu[i];
        for (mwSize k=0;k<m;++k) value+=S[i+m*k]*std::sqrt(std::max(eigenvalues[k],0.0))*noise[k];
        output[t+T*i]=value;
    }
}
void evaluate(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    require(nrhs==8&&nlhs==1,"Expected eight inputs and one output.");
    const double* Y=values(prhs[0]); const double* h=values(prhs[1]);
    const double* P=values(prhs[2]); const double* V=values(prhs[4]);
    const double* prior_mean=values(prhs[5]); const double* prior_cov=values(prhs[6]);
    const double* noise=values(prhs[7]);
    const mwSize T=mxGetM(prhs[0]),m=mxGetN(prhs[0]);
    require(T>=1&&m>=1&&m<=std::numeric_limits<mwSize>::max()/m,"Invalid panel dimensions.");
    const mwSize mm=m*m;
    require(mm<=std::numeric_limits<mwSize>::max()/T,"Panel exceeds the native indexing limit.");
    require(mxGetNumberOfDimensions(prhs[0])==2&&mxGetNumberOfDimensions(prhs[1])==2&&
        mxGetM(prhs[1])==T&&mxGetN(prhs[1])==m,"h must have the same shape as Y.");
    const mwSize* dimensions=mxGetDimensions(prhs[2]);
    require(dimensions[0]==m&&dimensions[1]==m&&mxGetNumberOfElements(prhs[2])==mm*T,
        "P must be m-by-m-by-T.");
    require(mxIsLogical(prhs[3])&&!mxIsSparse(prhs[3])&&
        mxGetNumberOfDimensions(prhs[3])==2&&mxGetM(prhs[3])==T&&mxGetN(prhs[3])==m,
        "mask must be a logical array with the same shape as Y.");
    require(mxGetM(prhs[4])==m&&mxGetN(prhs[4])==m&&mxGetNumberOfElements(prhs[4])==mm&&
        mxGetM(prhs[6])==m&&mxGetN(prhs[6])==m&&mxGetNumberOfElements(prhs[6])==mm,
        "V and prior covariance must be m-by-m.");
    require(mxGetNumberOfElements(prhs[5])==m&&mxGetM(prhs[7])==m&&mxGetN(prhs[7])==T&&
        mxGetNumberOfDimensions(prhs[7])==2,"Invalid prior mean or noise shape.");
    const mxLogical* mask=mxGetLogicals(prhs[3]);
    for (mwSize k=0;k<T*m;++k) {
        require(!mask[k]||(std::isfinite(Y[k])&&std::isfinite(h[k])),"Observed inputs must be finite.");
        require(std::isfinite(noise[k]),"Noise must be finite.");
    }
    for (mwSize k=0;k<mm*T;++k) require(std::isfinite(P[k]),"Correlations must be finite.");
    for (mwSize k=0;k<mm;++k) require(std::isfinite(V[k])&&std::isfinite(prior_cov[k]),"Covariances must be finite.");
    for (mwSize k=0;k<m;++k) require(std::isfinite(prior_mean[k]),"Prior mean must be finite.");
    std::vector<double> means(T*m),covariances(T*mm),mu(prior_mean,prior_mean+m),
        S(prior_cov,prior_cov+mm),predicted(mm),observed_covariance,observed_gain,
        innovation(m),conditional_mean(m),conditional_covariance(mm);
    std::vector<mwSize> observed;
    plhs[0]=mxCreateDoubleMatrix(T,m,mxREAL); double* B=mxGetDoubles(plhs[0]);
    for (mwSize t=0;t<T;++t) {
        predicted=S; if (t!=0) for (mwSize k=0;k<mm;++k) predicted[k]+=V[k];
        observed.clear(); for (mwSize i=0;i<m;++i) if (mask[t+T*i]) observed.push_back(i);
        const mwSize n=observed.size();
        if (n!=0) {
            observed_covariance.resize(n*n); observed_gain.resize(m*n);
            for (mwSize j=0;j<n;++j) {
                const mwSize oj=observed[j];
                innovation[j]=Y[t+T*oj]-mu[oj];
                for (mwSize i=0;i<m;++i) observed_gain[i+m*j]=predicted[i+m*oj];
                for (mwSize i=0;i<n;++i) {
                    const mwSize oi=observed[i];
                    observed_covariance[i+n*j]=predicted[oi+m*oj]+
                        std::exp(h[t+T*oi]/2)*std::exp(h[t+T*oj]/2)*P[mm*t+oi+m*oj];
                }
            }
            const auto gain=right_solve(observed_gain,observed_covariance,m,n);
            for (mwSize i=0;i<m;++i) for (mwSize j=0;j<n;++j) mu[i]+=gain[i+m*j]*innovation[j];
            for (mwSize j=0;j<m;++j) for (mwSize i=0;i<m;++i) {
                double value=predicted[i+m*j];
                for (mwSize k=0;k<n;++k) value-=gain[i+m*k]*predicted[observed[k]+m*j];
                S[i+m*j]=value;
            }
            symmetric(S,m);
        } else S=predicted;
        for (mwSize i=0;i<m;++i) means[t*m+i]=mu[i];
        std::copy(S.begin(),S.end(),covariances.begin()+mm*t);
    }
    normal_draw(S,mu,noise,B,T-1,T,m);
    for (mwSize step=1;step<T;++step) {
        const mwSize t=T-1-step;
        std::copy(covariances.begin()+mm*t,covariances.begin()+mm*(t+1),S.begin());
        predicted=S; for (mwSize k=0;k<mm;++k) predicted[k]+=V[k];
        const auto gain=right_solve(S,predicted,m,m);
        for (mwSize i=0;i<m;++i) {
            conditional_mean[i]=means[t*m+i];
            for (mwSize k=0;k<m;++k)
                conditional_mean[i]+=gain[i+m*k]*(B[t+1+T*k]-means[t*m+k]);
        }
        for (mwSize j=0;j<m;++j) for (mwSize i=0;i<m;++i) {
            double value=S[i+m*j];
            for (mwSize k=0;k<m;++k) value-=gain[i+m*k]*S[k+m*j];
            conditional_covariance[i+m*j]=value;
        }
        normal_draw(conditional_covariance,conditional_mean,noise+m*step,B,t,T,m);
    }
}
}
void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    try { evaluate(nlhs,plhs,nrhs,prhs); }
    catch (const Failure& failure) { mexErrMsgIdAndTxt(failure.id,"%s",failure.what()); }
    catch (const std::exception& failure) { mexErrMsgIdAndTxt("dsc:MeanNumerical","%s",failure.what()); }
}
