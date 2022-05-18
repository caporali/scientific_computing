function x = gauss_seidel_GR(W,ex,alpha,tol)
% Gauss-Seidel for GeneRank 
% Gauss-Sidel algorithm applied to a GeneRank problem
%
% 	input:
% 		W := GeneRank matrix 
% 		ex := constant term of the system
% 		alpha := damping factor
% 		tol := tolerance factor
% 	output
% 		x := estimate solution 
% 
% 	example
%		x = gauss_seidel_GR(W,ex,0.85,1e-14);
%
% Source: 
% Appunti del corso di Calcolo Scientifico
% Dario Andrea Bini e Lidia Aceto
% (important changes have been made)

fprintf("--- \n5. gauss_seidel_GR \n");

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

% computations: L, U
L = alpha*tril(W,-1);
U = speye(n) - triu(W)*diag(sparse(alpha*d_inv));

% short definitions
ex_complete = (1 - alpha)*ex;

% loop
tic;
while (res > tol)
	k = k + 2;

	x_old = x;

    x = d_inv.*x;
    x = L*x + ex_complete;
    % in case of non-symmetric matrix
    % x = L*x + alpha*sum(dang.*x) + ex_complete;
	x = U\x;
	
	res = norm(x - x_old,1);
	err = [err res res];
end

x = x/norm(x,1);

t = toc;

plot_function(alpha,5,err,k,t,x);