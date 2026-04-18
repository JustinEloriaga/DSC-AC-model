function [X, Y] = make_varXY(temp_YYY,nlag,intercept)

% temp_YYY  = nobs by nv matrix

% intercept = 1: include intercept
if nargin == 2
    intercept = 1;
end

% nlag = varinfo.nlag;
% nv   = varinfo.nv;

Y = [];
X = [];
for t = nlag+1:1:size(temp_YYY,1) %conditioning on first nlag obs
    
    Y = [Y; temp_YYY(t,:)];
    
    temp_XXX = [];
    for i = 1:1:nlag
        temp_XXX = [temp_XXX, temp_YYY(t-i,:)];
    end
    
    if intercept == 1
        temp_XXX = [temp_XXX, 1]; % add constant
    end
    
    X = [X; temp_XXX];
end