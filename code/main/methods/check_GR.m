function xreal = check_GR(W,ex,alpha)
% Check for GeneRank
% Direct check of GeneRank solutions 
%
% 	input:
% 		W := GeneRank matrix 
% 		ex := constant term of the system
% 		alpha := damping factor
%
%	output:
%		x_real := solution of GeneRank problem
% 
% 	example
%		xreal = check_GR(W,ex,alpha);

fprintf("--- \n0. check_GR \n");

ex = abs(ex);
ex = ex/norm(ex, 1); 

e = ones(size(W,1),1);
e = sparse(e);

d = sparse(W*e);
D = diag(d);
D_inv = diag(1./d);

% computation of M (complete matrix of GeneRank system)
tic;
M = speye(size(W,1)) - alpha*W*D_inv;
M = sparse(M);

% computation of xreal
xreal = M\((1 - alpha)*ex);
xreal = full(xreal/norm(xreal,1));

% text output: time and xreal
fprintf("\n	t = %f s \n", toc);
fprintf("	xreal(1:3,1) = [%e %e %e ...]' \n", xreal(1), xreal(2), xreal(3));

% check with  linear system (after resolution of dangling nodes)
A = alpha*W'*D_inv + (1 - alpha)*ex*e';
Axreal = A*xreal;
fprintf("	Axreal(1:3,1) = [%e %e %e ...]' \n", Axreal(1), Axreal(2), Axreal(3));
norm_diff(Axreal, xreal);