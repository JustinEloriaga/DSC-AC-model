function loglik = loglike_yt_given_ht_Pt(yt,ht,Pt,hhind,ht_prop)
% for yt ~ N(0,exp(ht/2)*Pt*exp(ht/2))


% update hhind-th columen of ht with ht_prop
ht(:,hhind) = ht_prop;

% size info
[m,~,T] = size(Pt);

% log-lik
loglik = 0;
% loglik_t = zeros(T,1);
MM = zeros(1,m);
for t=1:T
    VV = diag(exp(ht(t,:)/2))*Pt(:,:,t)*diag(exp(ht(t,:)/2));
    VV = (VV + VV')/2;
    temp_loglik = log_mvnpdf(yt(t,:), MM, VV);
    loglik = loglik + temp_loglik;
    
%     loglik_t(t) = temp_loglik;
end
