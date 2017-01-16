function [x,nref,lbl] = gen_clua(n,p,nc,sep)

%
% [x,n_obs_cl,lbl_cl] = gen_clua(num_obs,dimension,num_clusters,separation)
%
% Generation of clusters with normal observations and outliers
%
% Covariance matrices reasonably well conditioned
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

% Parameters

d_outl = 2;

xi299 = [ 6.635 9.210 11.345 13.277 15.086 16.812 18.475 20.090 21.666 23.209 ];
xi299 = [ xi299 24.725 26.217 27.688 29.141 30.578 32 33.409 34.805 36.191 37.566 ];
if p <= 20,
  vchi = xi299(p-1);
else
  vchi = (2.33 + sqrt(2*p-1))^2/2;
end
vchi = d_outl*vchi;

%% Scaling factors

fact1 = 5;        % largest scale for covariance matrix

% Number of observations in each cluster

nrf = floor(n*0.9);

nn = rand(nc,1);
n1 = sum(nn);
nref = floor(nn*(nrf-(p+1)*nc)/n1);
nref = nref + (p+1)*ones(nc,1);
n2 = sum(nref);
if n2 < nrf,
  nref(nc) = nref(nc) + nrf - n2;
end
nref = -sort(-nref);

% Labels for the observations in each cluster

lbl = zeros(nrf,1);
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
scl = fact1*sqrt(p);
for i = 1:nc,
  d = cntr(i,:);
  A = rand(p,p); A = A'*A;
  [U,D] = eig(A);
  v1 = scl*(rand(p,1) + 1.0e-3*ones(p,1));
  A = U*diag(sqrt(v1))*U';
  k = nref(i);
  x1 = normaliz(randn(k,p));
  x2 = x1*A;
  xaux = ones(k,1)*d + x2;
  x = [ x ; xaux ];
end

% Generation of outliers

k = n - nrf;
lc = 1;
y = x(find(lbl == lc),:);
m = mean(y); S = cov(y);
pass = 1;
while k > 0,
  if pass == 1,
    nat = floor(nref(lc)*0.1);
    if nat > k,
      nat = k;
    elseif nat < 1,
      nat = 1;
    end
    d = randn(1,p);
    uu = d*inv(S)*d';
    d = sqrt(vchi/uu)*d;
    xx = 0.1*randn(nat,p);
    xx = ones(nat,1)*(m + d) + xx;
    lbl = [ lbl ; -lc*ones(nat,1) ];
    x = [x ; xx ];
    k = k - nat;
    pass = 2;
  else
    d = randn(1,p);
    uu = d*inv(S)*d';
    d = sqrt(vchi/uu)*d;
    xx = m + d;
    lbl = [ lbl ; -lc ];
    x = [x ; xx ];
    k = k - 1;
    lc = lc + 1;
    if lc > nc,
      lc = 1;
    end
    y = x(find(lbl == lc),:);
    m = mean(y); S = cov(y);
    pass = 1;
  end
end
