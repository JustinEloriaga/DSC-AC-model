function [z_old, lik_old, n_try] = slice_sampling_v02(slice, z_old, lik_old)
% ---
% Elliptical slice sampling
% Generate z from p(z|data) where
% p(z|data) \propto L(data|z)*p(z)
% where 
% L(data|z): likelihood function
% p(z) = N(0, corr_z) 

% --------
% NOTE:
% _v02: p(z) = N(mu, cov_z) 
% slice.mean = mean(z)
% --------

% input: slice, z_old, lik_old
% slice.fcn_lik: function, likelihood function (in log)
% slice.nobs : integer, number of samples
% slice.chol_cov_z, nobs by nobs matrix, chol(cov_z)
% slice.scale_z: sometimes it is efficient to carry scale: chol( cov(z) ) = scale_z * chol_cov_z
% otherwise we set scale_z = 1 and use the full chol(cov(z))

% output: z_old, lik_old, n_try


% 1 Ellipse
% v = slice.chol_cov_z * randn(slice.nobs,1);

mu = slice.mean;

v = slice.chol_cov_z * randn(slice.nobs,1);
v = slice.scale_z*v;
% note v is de-meaned version of original v in the paper

z_old = z_old-mu; 
% note again that z_old is the de-meaned version of z_old;

% 2 Log-likelihood threshold
u = rand;
log_y = lik_old + log(u);

% 3 Initial proposal
theta = 0 + 2*pi*rand;

theta_min = theta-2*pi;
theta_max = theta;

% 4: Rest 
% z_new = (z_old-mu)*cos(theta) + (v-mu)*sin(theta) + mu;
z_new = z_old*cos(theta) + v*sin(theta) + mu;
lik_new = slice.fcn_lik(z_new);
valid = 0;
n_try = 1;

while (~valid)
    
    if lik_new > log_y
        % accept
        lik_old = lik_new;
        z_old = z_new;
        valid = 1;
    else
        %shrink the bracket and try a new point
        %disp('shrink ...');
        n_try = n_try + 1;
        
%         if n_try > 1000
%             keyboard
%         end
        
        if theta < 0
            theta_min = theta;
        else
            theta_max = theta;
        end
        theta = theta_min + (theta_max-theta_min)*rand;
        %z_new = (z_old-mu)*cos(theta) + (v-mu)*sin(theta) + mu;
        z_new = z_old*cos(theta) + v*sin(theta) + mu;
        lik_new = slice.fcn_lik(z_new);       
    end
	
end


