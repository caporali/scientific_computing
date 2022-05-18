function [v,H,m,mv] = arnoldi_process(W,deg,ex,v,d,m,mv)
% side-function of modified_arnoldi_GR.m 

w = ((1 - d)*sum(v))*ex + W*(d*(v./deg));

mv = mv + 1;
normw = norm(w);
H(1,1) = v'*w;
w = w - H(1, 1)*v;
h = norm(w);
if h < 0.618*normw
	h = v'*w;
	w = w - h*v;
	H(1,1) = H(1,1) + h;
	h = norm(w);
end
if h < 1e-14
	m = 1;
	H(2,1) = h;
	return;
else
	H(2,1) = h;
	w = w/h;
	v = [v w];
end
for j = 2:m
	w = ((1 - d)*sum(v(:,j)))*ex + W*(d*(v(:,j)./deg));
	
	mv = mv+1;
	normw = norm(w);
	H = [H;zeros(1,j-1)];
	if j < m
		H = [H zeros(j+1,1)];
	end
	for i = 1:j
		H(i,j) = v(:,i)'*w;
		w = w - H(i,j)*v(:,i);
	end
	h = norm(w);
	if h < 0.618 * normw 
		for i = 1:j
			h = v(:,i)'*w;
			w = w - h*v(:,i);
			H(i,j) = h + H(i,j);
		end
		h = norm(w);
	end
	if h < 1e-14
		m = j;
		H(j+1,j) = h;
		return;
	else
		H(j+1,j) = h;
		w = w/h;
		v = [v w];
	end
end

clear w;
clear h;