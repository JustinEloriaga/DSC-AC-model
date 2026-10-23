// Exact constant Gaussian mean update, with MATLAB-supplied standard normals.
// [mu,postmean,postcov] = dsc_constant_mean_mex(Y,h,P,mask,b0,S0,z)
// Only m-by-m workspaces are retained; missing dates contribute no information.
#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
struct Failure : std::runtime_error {
    const char* id;
    Failure(const char* identifier,const char* message)
        :std::runtime_error(message),id(identifier) {}
};
void require(bool condition,const char* id,const char* message) {
    if (!condition) throw Failure(id,message);
}
const double* values(const mxArray* a,const char* id) {
    require(mxIsDouble(a)&&!mxIsComplex(a)&&!mxIsSparse(a),id,
        "Numeric inputs must be real full doubles.");
    return mxGetDoubles(a);
}
bool vector_size(const mxArray* a,mwSize m) {
    return mxGetNumberOfDimensions(a)==2&&mxGetNumberOfElements(a)==m&&
        (mxGetM(a)==1||mxGetN(a)==1);
}
// Use a relative Frobenius symmetry test before averaging roundoff differences.
void validate_symmetry(const double* a,mwSize n,const char* id) {
    double norm=0.0,difference=0.0;
    for (mwSize j=0;j<n;++j) for (mwSize i=0;i<n;++i) {
        const double x=a[i+n*j],y=a[j+n*i];
        require(std::isfinite(x),id,"Covariance inputs must be finite.");
        norm=std::hypot(norm,x); difference=std::hypot(difference,x-y);
    }
    require(std::isfinite(difference)&&difference<=1e-12*std::max(1.0,norm),id,
        "Covariance inputs must be symmetric.");
}
void symmetrize(std::vector<double>& a,mwSize n) {
    for (mwSize j=0;j<n;++j) for (mwSize i=j+1;i<n;++i) {
        const double x=.5*a[i+n*j]+.5*a[j+n*i];
        a[i+n*j]=x; a[j+n*i]=x;
    }
}
void evaluate(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    const char* shape="dsc:ConstantMeanShape";
    const char* prior="dsc:ConstantMeanPrior";
    const char* covariance="dsc:ConstantMeanCovariance";
    const char* precision="dsc:ConstantMeanPrecision";
    const char* draw="dsc:ConstantMeanDraw";
    require(nrhs==7&&nlhs<=3,shape,"Expected seven inputs and at most three outputs.");
    const double* Y=values(prhs[0],shape);
    const double* h=values(prhs[1],shape);
    const double* P=values(prhs[2],shape);
    const mwSize T=mxGetM(prhs[0]),m=mxGetN(prhs[0]);
    require(mxGetNumberOfDimensions(prhs[0])==2&&T>=1&&m>=1&&
        m<=static_cast<mwSize>(std::numeric_limits<ptrdiff_t>::max())&&
        m<=std::numeric_limits<mwSize>::max()/m,shape,"Invalid panel dimensions.");
    const mwSize mm=m*m;
    require(T<=std::numeric_limits<mwSize>::max()/mm&&
        m<std::numeric_limits<mwSize>::max()&&
        m<=std::numeric_limits<mwSize>::max()/(m+1),shape,"Panel dimensions exceed indexing limits.");
    require(mxGetNumberOfDimensions(prhs[1])==2&&mxGetM(prhs[1])==T&&
        mxGetN(prhs[1])==m,shape,"h must have the same shape as Y.");
    const mwSize* pdims=mxGetDimensions(prhs[2]);
    require(mxGetNumberOfDimensions(prhs[2])<=3&&pdims[0]==m&&pdims[1]==m&&
        mxGetNumberOfElements(prhs[2])==mm*T,shape,"P must be m-by-m-by-T.");
    require(mxIsLogical(prhs[3])&&!mxIsSparse(prhs[3])&&
        mxGetNumberOfDimensions(prhs[3])==2&&mxGetM(prhs[3])==T&&
        mxGetN(prhs[3])==m,shape,"mask must be logical with the same shape as Y.");
    const mxLogical* mask=mxGetLogicals(prhs[3]);
    const double* b0=values(prhs[4],prior);
    const double* S0=values(prhs[5],prior);
    const double* z=values(prhs[6],draw);
    require(vector_size(prhs[4],m),prior,"The prior mean must be an m-vector.");
    require(mxGetNumberOfDimensions(prhs[5])==2&&mxGetM(prhs[5])==m&&
        mxGetN(prhs[5])==m,prior,"The prior covariance must be m-by-m.");
    require(vector_size(prhs[6],m),draw,"The normal draw must be an m-vector.");
    for (mwSize i=0;i<m;++i) {
        require(std::isfinite(b0[i]),prior,"The prior mean must be finite.");
        require(std::isfinite(z[i]),draw,"The normal draw must be finite.");
    }
    validate_symmetry(S0,m,prior);
    std::vector<double> factor(S0,S0+mm),Q(mm,0.0),rhs(m,0.0),postmean(m),
        observed_factor(mm),whitened(m*(m+1)),inverse_sd(m);
    std::vector<mwSize> observed(m);
    symmetrize(factor,m);
    const ptrdiff_t dimension=static_cast<ptrdiff_t>(m),one_column=1;
    const char lower='L',upper='U',left='L',normal='N',nonunit='N';
    const double one=1.0;
    ptrdiff_t info=0;
    dpotrf(&lower,&dimension,factor.data(),&dimension,&info);
    require(info==0,prior,"The prior covariance must be positive definite.");
    for (mwSize i=0;i<m;++i) Q[i+m*i]=1.0;
    dpotrs(&lower,&dimension,&dimension,factor.data(),&dimension,Q.data(),&dimension,&info);
    require(info==0,prior,"The prior precision solve failed.");
    for (mwSize i=0;i<m;++i) for (mwSize j=0;j<m;++j) rhs[i]+=Q[i+m*j]*b0[j];
    for (mwSize t=0;t<T;++t) {
        mwSize n=0;
        for (mwSize i=0;i<m;++i) if (mask[t+T*i]) observed[n++]=i;
        if (n==0) continue;
        for (mwSize j=0;j<n;++j) {
            const mwSize oj=observed[j],input=t+T*oj;
            require(std::isfinite(Y[input])&&std::isfinite(h[input]),
                "dsc:ConstantMeanObserved","Observed observations and log variances must be finite.");
            inverse_sd[j]=std::exp(-.5*h[input]);
            require(std::isfinite(inverse_sd[j])&&inverse_sd[j]>0,covariance,
                "An observed variance scale is outside the finite positive range.");
            for (mwSize i=0;i<n;++i)
                observed_factor[i+n*j]=P[mm*t+observed[i]+m*oj];
        }
        validate_symmetry(observed_factor.data(),n,covariance);
        for (mwSize j=0;j<n;++j) for (mwSize i=j+1;i<n;++i) {
            const double v=.5*observed_factor[i+n*j]+.5*observed_factor[j+n*i];
            observed_factor[i+n*j]=v; observed_factor[j+n*i]=v;
        }
        const ptrdiff_t no=static_cast<ptrdiff_t>(n),columns=no+1;
        dpotrf(&lower,&no,observed_factor.data(),&no,&info);
        require(info==0,covariance,"An observed correlation submatrix must be positive definite.");
        // [A,b] = chol(P_obs,'lower') \ [D_obs^-1,D_obs^-1*y_obs].
        // P alone is factored, so no time-varying covariance matrix is built.
        std::fill(whitened.begin(),whitened.begin()+n*(n+1),0.0);
        for (mwSize j=0;j<n;++j) {
            whitened[j+n*j]=inverse_sd[j];
            whitened[j+n*n]=inverse_sd[j]*Y[t+T*observed[j]];
        }
        dtrsm(&left,&lower,&normal,&nonunit,&no,&columns,&one,
            observed_factor.data(),&no,whitened.data(),&no);
        for (mwSize j=0;j<n;++j) {
            double r=0;
            for (mwSize k=0;k<n;++k) r+=whitened[k+n*j]*whitened[k+n*n];
            rhs[observed[j]]+=r;
            for (mwSize i=0;i<=j;++i) {
                double q=0;
                for (mwSize k=0;k<n;++k) q+=whitened[k+n*i]*whitened[k+n*j];
                Q[observed[i]+m*observed[j]]+=q;
            }
        }
    }
    // Only the upper triangle was accumulated; mirror it before validation.
    for (mwSize j=0;j<m;++j) for (mwSize i=j+1;i<m;++i) Q[i+m*j]=Q[j+m*i];
    for (double q:Q) require(std::isfinite(q),precision,"The posterior precision must be finite.");
    for (double r:rhs) require(std::isfinite(r),precision,"The posterior right-hand side must be finite.");
    dpotrf(&upper,&dimension,Q.data(),&dimension,&info);
    require(info==0,precision,"The posterior precision must be positive definite.");
    postmean=rhs;
    dpotrs(&upper,&dimension,&one_column,Q.data(),&dimension,postmean.data(),&dimension,&info);
    require(info==0,precision,"The posterior mean solve failed.");
    std::vector<double> noise(z,z+m);
    // Q=U'*U: U^-1*z has covariance U^-1*U'^-1 = Q^-1.
    dtrsm(&left,&upper,&normal,&nonunit,&dimension,&one_column,&one,
        Q.data(),&dimension,noise.data(),&dimension);
    for (mwSize i=0;i<m;++i)
        require(std::isfinite(postmean[i])&&std::isfinite(postmean[i]+noise[i]),draw,
            "The sampled mean must be finite.");
    if (nlhs>0) {
        plhs[0]=mxCreateDoubleMatrix(m,1,mxREAL);
        double* output=mxGetDoubles(plhs[0]);
        for (mwSize i=0;i<m;++i) output[i]=postmean[i]+noise[i];
    }
    if (nlhs>1) {
        plhs[1]=mxCreateDoubleMatrix(m,1,mxREAL);
        std::copy(postmean.begin(),postmean.end(),mxGetDoubles(plhs[1]));
    }
    if (nlhs>2) {
        plhs[2]=mxCreateDoubleMatrix(m,m,mxREAL);
        double* output=mxGetDoubles(plhs[2]);
        for (mwSize i=0;i<m;++i) output[i+m*i]=1.0;
        dpotrs(&upper,&dimension,&dimension,Q.data(),&dimension,output,&dimension,&info);
        require(info==0,precision,"The posterior covariance solve failed.");
        for (mwSize j=0;j<m;++j) for (mwSize i=j+1;i<m;++i) {
            const double v=.5*output[i+m*j]+.5*output[j+m*i];
            output[i+m*j]=v; output[j+m*i]=v;
        }
    }
}
}
void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    try { evaluate(nlhs,plhs,nrhs,prhs); }
    catch (const Failure& failure) { mexErrMsgIdAndTxt(failure.id,"%s",failure.what()); }
    catch (const std::exception& failure) {
        mexErrMsgIdAndTxt("dsc:ConstantMeanNumerical","%s",failure.what());
    }
}
