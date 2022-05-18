function x = modified_arnoldi_GR(W,ex,d,m,tol,maxmv)
% Modified Arnoldi method for GeneRank 
% Modified Arnoldi algorithm applied to a GeneRank problem
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
%		x = modified_arnoldi_GR(W,ex,0.85,1e-14);
%
% Source: 
% "Krylov subspace algorithms for computing GeneRank for the analysis of microarray data mining"
% Gang Wu, Wei Xu, Ying Zhang, Yimin Wei
% (no major changes)

fprintf("--- \n4. modified_arnoldi_GR \n");

ex = abs(ex);
ex = ex/sum(ex);

deg = sum(W);
ind = find(deg == 0);
deg(ind) = 1;
deg = deg';

% inizializations
w = ex;
err = []; 
mv = 0; 
r = 1;

% loop
tic;
while r > tol && mv <= maxmv
	mv1 = mv; 

	[v,H,m,mv] = arnoldi_process(W,deg,ex,w,d,m,mv); 
	if H(m+1,m) < 1e-14 
		[evr,ev] = eig(H(1:m,1:m));
		[~,ind] = max(diag(ev));
		v = v(:,1:m)*evr(:,ind); 
		return;
	end
	[u,s,t] = svd(H - [eye(m);zeros(1,m)],0);
	r = min(diag(s));  
	w = (((1 - d)*sum(v(:,m+1)))*ex + W*(d*(v(:,m+1)./deg))) - v(:,m+1);
	mv = mv+1;
	beta = w'*(v*u(:,m));
	gamma = norm(w)^2;
	beta = -r*beta/gamma;
	w = w + v(:,m+1);
	w = v*(H*t(:,m)) + beta*w;
	v = v*([t(:,m);beta]); 
	r = sum(abs(w-v))/sum(abs(v));

	err = [err r*ones(1, mv - mv1)];
end 

v = v/sum(abs(v));
v = abs(v);
err = err';

x = v;
k = mv;

t = toc;

plot_function(d,4,err,k,t,x);