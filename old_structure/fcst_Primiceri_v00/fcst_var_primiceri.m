function rst = fcst_var_primiceri(r, hmax, Yest, Yact, Xcalest, info)

% Similar to our "fcst_var_tvpsv.m"
% output of this function is exactly the same as "fcst_var_tvpsv.m"

% Note: For Primiceri estimation, Yest contains y including initial lags
% ----
% Primiceri and Del Negro's code (Something Minchul is not sure about)
% - I don't know why we put "-x" rather than "+x" in "tria(x)"
% ----
% Jonas's comment: I think this is because the definition of Z on page 845 in
% Primiceri (2005) equals -Z in the code. Hence, the code draws -aalpha



% Collect
y = Yest;

hmax = min(hmax, size(Yact,1));

lags = info.lags;
nburns = info.nburns;

r.B = r.B(:,:,nburns+1:end);
r.V = r.V(:,:,nburns+1:end);
r.A = r.A(:,:,nburns+1:end);
r.H = r.H(:,:,nburns+1:end);

%--------------- 

% Dimension
[T,k,ndraws] = size(r.B);
n = size(y,2);
na = size(r.A,2);
nh = size(r.H,2);

% Initialization
X1fix = [1];
for i=0:(lags-1)
    X1fix = [X1fix, y(end-i,:)];
end

% ---------------
% Step 1: Generate draws from predictive distribution

% matrix to store
Yfcst    = nan(hmax,n,ndraws);
Yfcst_st = nan(hmax,n,ndraws);

% predictive scores 
gb_plike_rao = nan(ndraws,hmax); %joint
gb_plike_rao_11 = nan(ndraws,hmax); % pred score for 1st variable
gb_plike_rao_22 = nan(ndraws,hmax); % pred score for 2nd variable
gb_plike_rao_33 = nan(ndraws,hmax); % pred score for 3rd variable
gb_plike_rao_44 = nan(ndraws,hmax); % pred score for 3rd variable


% iteration over simulation draws
record = 0;

for i_d = 1:1:ndraws
    
    % Parameter at i-th draw
    B0 = r.B(end,:,i_d);
    V0 = r.V(:,:,(i_d+1)); %+1 because the first is a prior draw
    A0 = r.A(end,:,(i_d+1)); %+1 because the first is a prior draw
    H0 = r.H(end,:,(i_d+1)); %+1 because the first is a prior draw

    % Variance of shocks to TVPs and SVs
    Vb0 = V0(1:k,1:k);
    Va0 = V0(k+1:k+na,k+1:k+na);
    Vh0 = V0(k+na+1:end,k+na+1:end);
    % correct loading: norm(blkdiag(Vb0,Va0,Vh0) - V0)
    
    sqVb0 = cholcov(Vb0);
    sqVa0 = cholcov(Va0);
    sqVh0 = cholcov(Vh0);
    
    X1 = X1fix;
    
    % B_stfcst = nan(m,n,hmax); %*do we have to do this?
    
    % Iteration for horizons
    for hind = 1:1:hmax
        
        % Move forward (Eqn (5-7) of Primiceri)
        H1 = H0 + randn(1,nh)*sqVh0;
        A1 = A0 + randn(1,na)*sqVa0;
        B1 = B0 + randn(1,k)*sqVb0;
        
        
        
        % Eqn (4) of Primiceri
        Y1 = kron(eye(n), X1)*B1' + ( tria(A1)\diag(exp(H1)) ) * randn(n,1);
        
        % Predictive score "conditional on B1,H1,A1
        fM = kron(eye(n), X1)*B1'; %mean
        fV = ( tria(A1)\diag(exp(H1)) ) * ( tria(A1)\diag(exp(H1)) )'; %variance
        
        gb_plike_rao(i_d,hind) = log(mvnpdf(Yact(hind,:), fM',fV));
        if gb_plike_rao(i_d,hind)==-Inf
            gb_plike_rao(i_d,hind) = -0.5*n*log(2*pi) - 0.5*LogAbsDet(fV) - 0.5*((Yact(hind,:)'-fM)'*(fV\eye(n))*(Yact(hind,:)'-fM));
        end
        gb_plike_rao_11(i_d,hind) = log(normpdf(Yact(hind,1), fM(1), sqrt(fV(1,1))));
        gb_plike_rao_22(i_d,hind) = log(normpdf(Yact(hind,2), fM(2), sqrt(fV(2,2))));
        gb_plike_rao_33(i_d,hind) = log(normpdf(Yact(hind,3), fM(3), sqrt(fV(3,3))));
        gb_plike_rao_44(i_d,hind) = log(normpdf(Yact(hind,4), fM(4), sqrt(fV(4,4))));

        % Update parameters
        H0 = H1;
        A0 = A1;
        B0 = B1;
        X1 = [1, Y1', X1(2:(end-n))]; % always include intercept
        
        % Store
        Yfcst(hind,:,i_d) = Y1;
        %B_stfcst(:,:,hind) = B1 %*if we want to do this, we have to reshape B1

    end
    
    % NON-STATIONARITY CHECK, DO WE WANT TO DO THIS FOR THIS MODEL AS WELL? 

end

% ----------------------------
% Step 2: Evaulate forecasts

% Mean fcst
Mfcst  = mean(Yfcst,3);
Mefcst = median(Yfcst,3);

e     = Yact - Mfcst; %point error

% Interval fcst (70%)
IfcstL70 = quantile(Yfcst,0.15,3);
IfcstU70 = quantile(Yfcst,0.85,3);
i70 = (IfcstL70 < Yact) & (IfcstU70 > Yact); %hit-indicator

% Interval fcst (80%)
IfcstL80 = quantile(Yfcst,0.1,3);
IfcstU80 = quantile(Yfcst,0.9,3);
i80 = (IfcstL80 < Yact) & (IfcstU80 > Yact); %hit-indicator

% Interval fcst (90%)
IfcstL90 = quantile(Yfcst,0.05,3);
IfcstU90 = quantile(Yfcst,0.95,3);
i90 = (IfcstL90 < Yact) & (IfcstU90 > Yact); %hit-indicator


% ----------------------------
% Step 3: Collect results
rst = [];

% forecasts
rst.act    = Yact;
rst.mean   = Mfcst;    % mean forecast 
rst.median = Mefcst;   % mean forecast 
rst.L70    = IfcstL70; % interval forecasts
rst.U70    = IfcstU70;
rst.L80    = IfcstL80;
rst.U80    = IfcstU80;
rst.L90    = IfcstL90;
rst.U90    = IfcstU90;

% error statistics
rst.e   = e;   %point error
rst.i70 = i70; %interval hit indicators
rst.i80 = i80;
rst.i90 = i90;


% other useful objects
rst.Xcalest = Xcalest;
rst.Yest    = Yest;

rst.plike_rao_mean      = mean(exp(gb_plike_rao),1);
rst.plike_rao_11_mean   = mean(exp(gb_plike_rao_11),1);
rst.plike_rao_22_mean   = mean(exp(gb_plike_rao_22),1);
rst.plike_rao_33_mean   = mean(exp(gb_plike_rao_33),1);
rst.plike_rao_44_mean   = mean(exp(gb_plike_rao_44),1);

rst.log_plike_rao_mean      = mean(gb_plike_rao,1);
rst.log_plike_rao_11_mean   = mean(gb_plike_rao_11,1);
rst.log_plike_rao_22_mean   = mean(gb_plike_rao_22,1);
rst.log_plike_rao_33_mean   = mean(gb_plike_rao_33,1);
rst.log_plike_rao_44_mean   = mean(gb_plike_rao_44,1);


