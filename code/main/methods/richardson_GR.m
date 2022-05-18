function x = richardson_GR(W,ex,alpha,tol)
% Richardson for GeneRank 
% Richardson algorithm applied to a GeneRank problem
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
%		x = richardson_GR(W,ex,0.85,1e-14);
%
% Source: 
% Appunti del corso di Calcolo Scientifico
% Dario Andrea Bini e Lidia Aceto
% (important changes have been made)

fprintf("--- \n7. richardson_GR \n");

% beta
beta = 0.7;

n = size(W,1);

ex = abs(ex);
ex = ex/norm(ex,1);   

e = ones(n,1);
d = W*e;
d_inv = 1./d;

% inizializations: k, x, res, err
k = 0;	
x = e/n;
res = 1;
err = [];

% short definitions
ex_complete = beta*(1 - alpha)*ex;
W_complete = beta*alpha*W;

% loop
tic;
while (res > tol)
	k = k + 1;

	x_old = x;

	x = (1 - beta)*x + W_complete*(d_inv.*x) + ex_complete;
	
	res = norm(x - x_old,1);
	err = [err res];
end

x = x/norm(x,1);

t = toc;

plot_function(alpha,7,err,k,t,x);