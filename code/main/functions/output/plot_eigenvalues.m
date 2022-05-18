function plot_eigenvalues(W,alpha)
% Plot eigenvalues of GeneRank system and preconditioned GeneRank system 
%
% 	input:
% 		W := GeneRank matrix 
% 		alpha := damping factor
% 
% 	example
%		plot_eigenvalues(W,0.85)

% remove dangling nodes
W = remove_dangling(W);

n = size(W,1);

e = ones(n,1);
d = W*e;
% in case of non-symmetric matrix
% dang = d==0;
% d = d + dang*n;
D = diag(d);
d_inv = 1./d;
D_inv = diag(d_inv);
rad_D_inv = sqrtm(D_inv);

M = D - alpha*W;
pM = eye(n) - alpha*rad_D_inv*W*rad_D_inv; 

eig_M = sort(abs(eig(M)), "descend");
eig_pM = sort(abs(eig(pM)), "descend");

figure;
plot(eig_M, "g.");
title("$$D - \alpha W$$","interpreter","latex");
xlabel("n","interpreter","latex");
ylabel("eigenvalue","interpreter","latex");
ax=gca;
ax.XAxis.FontSize = 8;
ax.YAxis.FontSize = 8;
ax.XLabel.FontSize = 10;
ax.YLabel.FontSize = 10;
ax.Title.FontSize = 12;

figure;
plot(eig_pM, "r.");
title("$$I - \alpha D^{\frac{-1}{2}} W D^{\frac{-1}{2}}$$","interpreter","latex");
xlabel("n", "interpreter", "latex");
ylabel("eigenvalue", "interpreter", "latex");
ax=gca;
ax.XAxis.FontSize = 8;
ax.YAxis.FontSize = 8;
ax.XLabel.FontSize = 10;
ax.YLabel.FontSize = 10;
ax.Title.FontSize = 12;