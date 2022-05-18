function x = pcg_GR(W,ex,alpha,tol)
% PCG for GeneRank 
% Preconditioned conjugate gradient algorithm applied to a GeneRank problem
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
%		x = pcg_GR(W,ex,0.85,1e-14);
%
% Source: 
% "A preconditioned conjugate gradient algorithm for GeneRank with application to microarray data mining"
% Gang Wu, Wei Xu, Ying Zhang, Yimin Wei
% (important changes have been made)

fprintf("--- \n2. pcg_GR \n");

n = size(W,1);

ex = abs(ex);
ex = ex/norm(ex, 1);   

e = ones(n,1);
d = W*e;
d_inv = 1./d;

% inizializations: k, x, r, p, y, res, err
k = 0;	
x = zeros(n, 1);
r = ex;
p = ex.*d_inv;
y = p;
res = 1;
err = [];

% short definitions
W_complete = alpha*W; 

% loop
tic;
while (res > tol)
	k = k + 1;

	y_old = y;
	r_old = r;
	
	% first loop, inizializations: z, nu, mu
	z = d.*p - W_complete*p;
	nu = (y_old'*r_old)/(p'*z);
	x = x + nu*p;
	r = r_old - nu*z;
	y = r.*d_inv;
	mu = (y'*r)/(y_old'*r_old);
	p = y + mu*p;

	res = norm(r,1);
	err = [err res];
end

x = d.*x;
x = x/norm(x,1);

t = toc;

plot_function(alpha,2,err,k,t,x);