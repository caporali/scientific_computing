function W = remove_dangling(W)
% Remove dangling nodes
%
% 	input:
% 		W := GeneRank matrix 
%
%   output:
% 		W := GeneRank matrix without dangling nodes
% 
% 	example
%		W = remove_dangling(W);
    
e = ones(size(W,1),1);
e = sparse(e);

d = W*e;
dang = d==0;
dang = sparse(dang);

Dang = dang*e';
Dang_t = Dang' - diag(diag(Dang'));

W = W + Dang + Dang_t;
W = sparse(W);