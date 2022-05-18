function A = random_GR_matrix(n,k)
% Random GeneRank matrix
% Generate a boolean symmetric matrix of size n and density k/n. Possible input of GeneRank problem
%
%	input: 
%		n := size of the matrix
%       k := density (k/n)
%
%	output:
%		A := matrix for GeneRank
%
%	example:
%		A = random_GR_matrix(10000,7);

A = sprandsym(n,k/n);
A = A ~= 0;
A = A - diag(diag(A));