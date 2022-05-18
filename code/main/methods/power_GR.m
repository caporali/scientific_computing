function x = power_GR(W,ex,alpha,tol)
% Power method for GeneRank 
% Power algorithm applied to a GeneRank problem
%
%	input:
% 		W := GeneRank matrix 
% 		ex := constant term of the system
% 		alpha := damping factor
% 		tol := tolerance factor
%	output
% 		x := estimate solution 
% 
%	example
%         x = power_GR(W,ex,0.85,1e-14);
% Source: 
% Appunti del corso di Calcolo Scientifico
% Dario Andrea Bini e Lidia Aceto
% (important changes have been made)

fprintf("--- \n6. power_GR \n");

n = size(W,1);

ex = abs(ex);
ex = ex/norm(ex, 1);   

e = ones(n,1);
d = W*e;
% in case of non-symmetric matrix
% dang = d==0;
% d = d + dang*n;
d_inv = 1./d;

% inizializations: k, x, res, err
k = 0;	
x = e/n;
res = 1;
err = [];

% short definitions
ex_complete = (1 - alpha)*ex;

% loop ( A = alpha*W*D_inv  + (1 - alpha)*ex*ones(n,1)' )
tic;
while res > tol
	k = k + 1;

	x_old = x;

    x = d_inv.*x;
    x = W*x;
    % in case of non-symmetric matrix
    % x = W*x + sum(dang.*x);
    x = alpha*x + ex_complete;

    res = norm(x - x_old, 1);
	err = [err res];
end

x = x/norm(x,1);

t = toc;

plot_function(alpha,6,err,k,t,x);