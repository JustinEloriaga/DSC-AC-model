function rst = fcst_var_primiceri_msv2(r, hmax, Yest, Yact, Xcalest, info)

% Similar to our "fcst_var_tvpsv.m"

% msv2: for Hansen's model

% Collect
y = Yest;

hmax = min(hmax, size(Yact,1));

lags   = info.lags;

nburns = info.nburns / info.nthin;

r.B = r.B(:,:,nburns+1:end);
r.V = r.V(:,:,nburns+1:end);

r.r = r.r(:,:,nburns+1:end);
r.sig2r = r.sig2r(nburns+1:end,:);
r.P = r.P(:,:,:,nburns+1:end);

r.h = r.h(:,:,nburns+1:end);
r.sig2h = r.sig2h(nburns+1:end,:);


%---------------

% Dimension
[T,k,nsave] = size(r.B);
ndraws = nsave * info.nthin;

n = size(y,2);
n_r = n*(n-1)/2; % # of correlation elements
% na = size(r.A,2);
% nh = size(r.H,2);

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

count = 1;
for i_d = 1:nsave
    
    for i_e = 1:info.nthin
        
        % Parameter at i-th draw
        B0 = r.B(end,:,i_d);
        V0 = r.V(:,:,(i_d)); %+1 because the first is a prior draw
        
        % Correlation
        %Q0 = r.Q(:,:,end,i_d);
        %d0 = r.d(i_d,1);
        %k0 = r.k(i_d,1);
        %A0 = eye(n) / r.Ainv(:,:,i_d);
        
        r0 = r.r(end,:,i_d);
        P0 = r.P(:,:,end,i_d);
        sig2r0 = r.sig2r(i_d,:);
        
        % Volatility
        H0 = r.h(end,:,i_d);
        sig2h0 = r.sig2h(i_d,:);
        
        
        % Variance of shocks to TVPs and SVs
        Vb0 = V0(1:k,1:k);
        sqVb0 = cholcov(Vb0);
        
        X1 = X1fix;
        
        % B_stfcst = nan(m,n,hmax); %*do we have to do this?
        
        % Iteration for horizons
        for hind = 1:1:hmax
            
            % Move forward (Eqn (5-7) of Primiceri)
            
            % time-varying parameter
            B1 = B0 + randn(1,k)*sqVb0;
            
            % volatility
            H1 = H0 + sqrt(sig2h0).*randn(1,n); %variance
            
            % correlation
            r1 = r0 + sqrt(sig2r0).*randn(1,n_r); %correlation
            
%             P1 = veclAtoC_4v_mex(r0);
            P1 = veclAtoC(r0);
            
            % Eqn (4) of Primiceri
            %Y1 = kron(eye(n), X1)*B1' + ( tria(A1)\diag(exp(H1)) ) * randn(n,1);
            Y1 = kron(eye(n), X1)*B1' + diag(exp(H1/2))*chol(P1,'lower') * randn(n,1);
            
            % Predictive score "conditional on B1,H1,A1
            fM = kron(eye(n), X1)*B1'; %mean
            %fV = ( tria(A1)\diag(exp(H1)) ) * ( tria(A1)\diag(exp(H1)) )'; %variance
            fV =  diag(exp(H1/2))*P1* diag(exp(H1/2));
            fV = (fV + fV')/2;
            
            gb_plike_rao(count,hind) = log(mvnpdf(Yact(hind,:), fM',fV));
            if gb_plike_rao(count,hind)==-Inf
                gb_plike_rao(count,hind) = -0.5*n*log(2*pi) - 0.5*LogAbsDet(fV) - 0.5*((Yact(hind,:)'-fM)'*(fV\eye(n))*(Yact(hind,:)'-fM));
            end
            gb_plike_rao_11(count,hind) = log(normpdf(Yact(hind,1), fM(1), sqrt(fV(1,1))));
            gb_plike_rao_22(count,hind) = log(normpdf(Yact(hind,2), fM(2), sqrt(fV(2,2))));
            gb_plike_rao_33(count,hind) = log(normpdf(Yact(hind,3), fM(3), sqrt(fV(3,3))));
            gb_plike_rao_44(count,hind) = log(normpdf(Yact(hind,4), fM(4), sqrt(fV(4,4))));
            
            % Update parameters
            r0 = r1;
            H0 = H1;
            B0 = B1;
            X1 = [1, Y1', X1(2:(end-n))]; % always include intercept
            
            % Store
            Yfcst(hind,:,count) = Y1;
            %B_stfcst(:,:,hind) = B1 %*if we want to do this, we have to reshape B1
            
        end
        
        % NON-STATIONARITY CHECK, DO WE WANT TO DO THIS FOR THIS MODEL AS WELL?
        
        count = count + 1;
    end
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

rst.sd  = std(Yfcst,[],3); % standard deviation
rst.iqr = iqr(Yfcst, 3); % IQR

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

% ---
% entire predictive density
if info.save_pred_dens == 1
    rst.Yfcst=Yfcst;
end

% ---
% CRPS
rst.crps = size(Yact);
for vvind = 1:size(Yact,1)
    for hhind = 1:size(Yact,2)
        rst.crps(vvind,hhind) = compute_crps2( Yact(vvind,hhind) , squeeze(Yfcst(vvind,hhind,:)) );
    end
end


