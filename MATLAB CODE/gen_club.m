function [x,nref,lbl,flg] = gen_club(n,p,nc,sep,mode)

%
% [x,n_obs_cl,lbl_cl] = gen_club(num_obs,dimension,num_clusters,separation,mode)
%
% Generation of clusters. Not too ill conditioned
% mode = 1,  normals
% mode = 2,  uniforms
% mode = 3,  student-t (p degrees of freedom)
%

% DP/FJP  6/29/01

if p > 100,
  disp('Dimension too high');
  return
end

if (nc*(p+1) > n)|(5*p > n),
  disp('Number of obseravtions too small');
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
%% Scale factors

fact1 = 5;        % Largest scale for the covariance matrix
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

% Labels identifying each cluster

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

% Generation of the centers of the clusters

scl = sep/sqrt(p);
cntr = [];
for i = 1:nc,
  d = randn(1,p);
  cntr = [ cntr ; (scl*d) ];
end

% Generation of the observations in each cluster

x = [];
scl = fact1*sqrt(p);
for i = 1:nc,
  d = cntr(i,:);
  A = rand(p,p); A = A'*A;
  [U,D] = eig(A);
  v1 = scl*(rand(p,1) + 1.0e-3*ones(p,1));
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