function [C , iter_number ] = veclAtoC_paper( gamma_in )
% take from the original paper


% ============
% tol_value = 1e-8 * sqrt(4); % as in the original paper
tol_value = 1e-8 ; % as in the original paper
% ============

C = [];
iter_number = -1;
% Check if input is of proper format : gamma is of suitable length
% and tolerance value belongs to a proper interval
n = 0.5*(1+ sqrt (1+8* length ( gamma_in )));
if isvector ( gamma_in ) && ( n == floor ( n ))...
        && ( tol_value > 0 && tol_value < 1)
    % Place elements from gamma into off - diagonal parts
    % and put zeros on the main diagonal of nxn symmetric matrix A
    A = zeros ( n );
    A ( logical ( tril ( ones ( n ) , -1))) = gamma_in ;
    A = A + A';
    % Read some properties of matrix A
    diag_vec = diag ( A );
    diag_ind = logical ( eye ( n ));
    % Iterative algorithm to get the proper diagonal vector
    dist = sqrt ( n );
    while dist > sqrt ( n )* tol_value
        diag_delta = log ( diag ( expm ( A )));
        diag_vec = diag_vec - diag_delta ;
        A ( diag_ind ) = diag_vec ;
        dist = norm ( diag_delta );
        iter_number = iter_number + 1;
    end
    % Get a unique reciprocal correlation matrix
    C = expm ( A );
    C ( diag_ind ) = ones (n ,1);
else
    fprintf ('Error : input is of wrong format');
end