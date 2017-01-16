function val = sim_clus(n,p,k,sep,num,mode,cutoff)

%
% val = sim_clus(num_obs,dimension,num_clusters,
%                  separation,num_replications,mode,cutoff)
%
% Simulation experiment for the kurtosis cluster identification
% procedure
%
% val : probability of failure [ < 50% / mix / both ]
%       mean number of failures / max eig cov within groups
%       mix errors only if more than 5% observations are mixed
%
% mode :  0 , normals with outliers
%        >0 , normals(1)/uniforms(2)/student-t(3) ill-conditioned
%        <0 , normals(-1)/uniforms(-2)/student-t(-3) better conditioned
%

% DP/FJP  6/29/01

if nargin < 7,
  cutoff = 0.1;
end
if nargin < 6,
  mode = 1;
end
if nargin < 5,
  num = 10;
end
if nargin < 4,
  sep = 4;
end

v5 = 0;
vm = 0;
va = 0;
nerr = 0;
rff1 = 0;
rff2 = 0;

for i = 1:num,

  if mode < 0,
    [x,nref0,lbl0] = gen_club(n,p,k,sep,abs(mode));
  elseif mode == 0,
    [x,nref0,lbl0] = gen_clua(n,p,k,sep);
  else
    [x,nref0,lbl0] = gen_clus(n,p,k,sep,abs(mode));
  end
  nref = [ 0 ; cumsum(nref0) ];

  lbl = clus_kur(x);
  k1 = max(lbl);

  if mode == 0,
    kmn = abs(min(lbl0));
    keff = k + kmn;
  else
    keff = k;
  end
  lrf = zeros(keff,k1);
  vv3 = lrf;
  rf1 = 0;
  rf2 = 0;
  for j = 1:k,
    lblfd = find(lbl0 == j);
    lbl1 = lbl(lblfd);
    vv3(j,:) = hist(lbl1,1:k1);
    lrf(j,:) = vv3(j,:)/sum(lbl0 == j);
  end
  if mode == 0,
    for j = 1:kmn,
      lblfd = find(lbl0 == -j);
      lbl1 = lbl(lblfd);
      vv3(j+k,:) = hist(lbl1,1:k1);
      lrf(j+k,:) = vv3(j+k,:)/sum(lbl0 == -j);
    end
    lrfs = zeros(k,k1);
    for j = 1:k,
      lblfd = find(abs(lbl0) == j);
      lbl1 = lbl(lblfd);
      lrfs(j,:) = hist(lbl1,1:k1);
    end
  end

  vv1 = max(lrf');
  if min(vv1) <= 0.5,
    rf1 = 1;
    v5 = v5 + 1;
  end
  vv2 = sum(lrf > 0.05);
  if (sum(vv2 > 1) > 0),
    rf2 = 1;
    vm = vm + 1;
  end
  if rf1 + rf2 > 1,
    va = va + 1;
  end
  if k + k1 < 3,
    vv3 = [ vv3 zeros(k,1) ];
  end
  if keff <= k1,
    [vv4,ix4] = max(vv3');
  else
    [vv4,ix4] = max(vv3);
  end
  ref1 = sum(vv4);
  if mode == 0,
    if k <= k1,
      [vv4,ix4] = max(lrfs');
    else
      [vv4,ix4] = max(lrfs);
    end
    ref1 = max(ref1,sum(vv4));    
  end
  nerr = nerr + n - ref1;

% Within clusters variability computation

  S = zeros(p,p);
  for j = 1:k1,
    aix = (lbl == j);
    idx = find(aix);
    nx = sum(aix);
    if nx > 1,
      xx = x(idx,:);
      S = S + (nx-1)*cov(xx);
    end
  end
  Sc = zeros(p,p);
  for j = 1:max(lbl0),
    aix = (lbl0 == j);
    idx = find(aix);
    nx = sum(aix);
    if nx > 1,
      xx = x(idx,:);
	Sc = Sc + (nx - 1)*cov(xx);
    end
  end
  S1 = (n-1)*cov(x);

  rff1 = rff1 + max(eig(S))/max(eig(Sc));
  rff2 = rff2 + sum(diag(S))/sum(diag(Sc));

end

val = [v5 vm va nerr rff1 rff2]/num;
