function x = jacobi_GR(W,ex,alpha,tol)
% Jacobi for GeneRank 
% Jacobi algorithm applied to a GeneRank problem
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
%		x = jacobi_GR(W,ex,0.85,1e-14);
%
% Source: 
% "A preconditioned conjugate gradient algorithm for GeneRank with application to microarray data mining"
% Gang Wu, Wei Xu, Ying Zhang, Yimin Wei
% (important changes have been made)

fprintf("--- \n3. jacobi_GR \n");

n = size(W,1);

ex = abs(ex);
ex = ex/norm(ex, 1);   

e = ones(n,1);
d = W*e;
d_inv = 1./d;

% inizializations: k, x, res, err
k = 0;	
x = e/n;
res = 1;
err = [];

% short definitions
ex_complete = (1 - alpha)*ex;
W_complete = alpha*W;

% loop
tic;
while (res > tol)
	k = k + 1;

	x_old = x;
	
	x = W_complete*(d_inv.*x) + ex_complete;
	
	res = norm(x - x_old,1);
	err = [err res];
end

x = x/norm(x,1);

t = toc;

plot_function(alpha,3,err,k,t,x);