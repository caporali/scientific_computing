function norm_diff(x,xreal)
% Print norm of residue
%
% 	input:
% 		x := estimate solution 
% 		xreal := solution of GeneRank problem
% 
% 	example
%		norm_diff(x,xreal);

fprintf("	||x - xreal||_1 = %e\n", norm(abs(x - xreal), 1));