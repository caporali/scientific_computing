function sin_angle(x, xreal)
% Compute and print $\sin\angle(x, xreal)$
%
% 	input:
% 		x := estimate solution 
% 		xreal := solution of GeneRank problem
% 
% 	example
%		sin_angle(x, xreal);

xnorm = x/norm(x,2);
xrealnorm = xreal/norm(xreal,2);

tmp = xnorm' * xrealnorm;
tmp = xnorm * tmp;

s = norm(xrealnorm - tmp,2);

fprintf("	sin angle(x, xreal) = %e\n", s);