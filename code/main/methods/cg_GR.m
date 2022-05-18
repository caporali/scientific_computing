function [x, err] = cg_GR(W,ex,alpha,tol)
% CG for GeneRank 
% Conjugate gradient algorithm applied to a GeneRank problem
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
%		x = cg_GR(W, ex, 0.85, 1e-14);
%
% Source: 
% "A preconditioned conjugate gradient algorithm for GeneRank with application to microarray data mining"
% Gang Wu, Wei Xu, Ying Zhang, Yimin Wei
% (important changes have been made)

fprintf("--- \n1. cg_GR \n");

n = size(W,1);

ex = abs(ex);
ex = ex/norm(ex, 1);   

e = ones(n,1);
d = W*e;

% inizializations: k, x, r, p, res, err
k = 0;	
x = zeros(n, 1);
r = ex;
p = ex;
res = 1;
err = [];

% short definitions
W_complete = alpha*W; 

% loop
tic;
while (res > tol)
	k = k + 1;

	r_old = r;
	
	% first loop, inizializations: z, nu, mu
	z = d.*p - W_complete*p;
	nu = (r_old'*r_old)/(p'*z);
	x = x + nu*p;
	r = r_old - nu*z;
	mu = (r'*r)/(r_old'*r_old);
	p = r + mu*p;

	res = norm(r,1);
	err = [err res];
end

x = d.*x;
x = x/norm(x,1);

t = toc;

plot_function(alpha,1,err,k,t,x);