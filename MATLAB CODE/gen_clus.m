function [x,nref,lbl,flg] = gen_clus(n,p,nc,sep,mode)

%
% [x,n_obs_cl,lbl_cl] = gen_clus(num_obs,dimension,num_clusters,separation,mode)
%
% Generation of clusters. Ill-conditioned covariance matrices
% mode = 1,  normals
%      = 2,  uniforms
%      = 3,  student-t (p degrees of freedom)
%

% DP/FJP  6/29/01

if p > 100,
  disp('Dimensión too high');
  return
end

if (nc*(p+1) > n)|(5*p > n),
  disp('Too few observations');
  return
end

if sep < 0.05,
  disp('Separation too small');
  return
end

if mode == 3,
  ngl = p;
end

% Parameters
%% Scaling factors

fact1 = 10;        % largest scale for covariance matrices
fctuni = (3/2)^(1/3);

% Number of observations in each cluster

nn = rand(nc,1);
n1 = sum(nn);
nref = floor(nn*(n-(p+1)*nc)/n1);
nref = nref + (p+1)*ones(nc,1);
n2 = sum(nref);
if n2 < n,
  nref(nc) = nref(nc) + n - n2;
end

% Labels for the observations in each cluster

lbl = zeros(n,1);
i = 1;
lbl1 = 1;
lbl2 = 0;
while (i <= nc),
  lbl2 = nref(i) + lbl2;
  lbl(lbl1:lbl2) = i*ones((lbl2-lbl1+1),1);
  lbl1 = lbl2 + 1;
  i = i + 1;
end

% Generation of the centers of the distributions

scl = sep/sqrt(p);
cntr = [];
for i = 1:nc,
  d = randn(1,p);
  cntr = [ cntr ; (scl*d) ];
end

% Generation of the data in each cluster

x = [];
scl = fact1*p^2;
scl1 = scl^(1/(p-1));
v1 = scl1.^(0:(p-1));
v1 = v1./sqrt(scl);
for i = 1:nc,
  d = cntr(i,:);
  A = rand(p,p); A = A'*A;
  [U,D] = eig(A);
  A = U*diag(sqrt(v1))*U';
  k = nref(i);
  if mode == 1,
    x1 = normaliz(randn(k,p));
    x2 = x1*A;
  elseif mode == 2,
    x1 = -fctuni + 2*fctuni*rand(k,p);
    x2 = x1*A;
  elseif mode == 3,
    x1 = normaliz(randn(k,p));
    x3 = x1*A;
    x4 = norm(randn(ngl,1))/sqrt(ngl);
    x2 = x3./x4;
  end
  xaux = ones(k,1)*d + x2;
  x = [ x ; xaux ];
end
