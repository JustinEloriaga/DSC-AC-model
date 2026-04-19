function loglik = loglike_yt_given_rt_ht(yt,ht,rt,rrind,rt_prop)
% for yt ~ N(0,exp(ht/2)*Pt*exp(ht/2))


% update hhind-th columen of ht with ht_prop
rt(:,rrind) = rt_prop;

% size info
% [m,~,T] = size(Pt);
[T,m] = size(ht);

% from rt to Pt(correlation matrix)
Pt = zeros(m,m,T);
for t=1:T
%     Pt(:,:,t) = veclAtoC_4v_mex(rt(t,:)); %4v works for 4 variables
    Pt(:,:,t) = veclAtoC(rt(t,:));
end

% log-lik
loglik = 0;
% loglik_t = zeros(T,1);
MM = zeros(1,m);
for t=1:T
    VV = diag(exp(ht(t,:)/2))*Pt(:,:,t)*diag(exp(ht(t,:)/2));
    VV = (VV + VV')/2;
    
    try
        temp_loglik = log_mvnpdf(yt(t,:), MM, VV);
        loglik = loglik + temp_loglik;
    catch
        disp('**** given_rt_ht, log_mvnpdf error ...');
        loglik = -Inf;
        return
    end
%     loglik_t(t) = temp_loglik;
end
