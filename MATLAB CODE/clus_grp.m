function [lbl1,ncl1] = clus_grp(x,lbl)

%
% [label_r,nclus_r] = clus_grp(x,label)
%
% Cluster identification from projections onto directions
% maximizing and minimizing the kurtosis coefficient of
% the data
%
% Routine to test for the merging of clusters (to be used as
% a subroutine of clus_kur)
%
% Inputs:   observations, x, matrix with one observation in each row
%           label,  initial labels
% Outputs:  label_r,  revised labels
%           nclus_r,  number of revised clusters
%

% DP/FJP  6/29/01

% Initial data

[n,p] = size(x);
[lbl0,ncl0] = ord_clus(lbl);

if ncl0 == 1,
  lbl1 = lbl0;
  return
end

% Percentiles for xi2 distribution

xi299 = [ 6.635 9.210 11.345 13.277 15.086 16.812 18.475 20.090 21.666 23.209 ];
xi299 = [ xi299 24.725 26.217 27.688 29.141 30.578 32 33.409 34.805 36.191 37.566 ];

lbl1 = lbl0;

% Ordering clusters by size

aa = [];
k = 1;
while k <= ncl0,
  aa = [ aa sum(lbl0 == k) ];
  k = k + 1;
end
[ab,ib] = sort(-aa);      % direct order
[dummy,iz] = sort(ib);    % inverse order

% Cutoff value

if p <= 20,
  xref = xi299(p-1);
else
  xref = (2.33 + sqrt(2*p-1))^2/2;
end

% Consider clusters from largest to smallest
% Merge if appropriate

k = 1;
while k <= ncl0,
  ic = ib(k);
  rv = (sum(lbl1 == ic) > p);
  while rv,
    rv = 0;
    xx = x(find(lbl1 == ic),:);
    m = mean(xx);
    S = cov(xx);
    S1 = inv(S);

    ix = find(lbl1 ~= ic);
    lix = length(ix);
    if lix > 0,
      xy = x(ix,:) - ones(lix,1)*m;
      xz = xy*S1;
      dd = sum((xz.*xy)')';
      iw = ix(find(dd <= xref));
      liw = length(iw);
      if liw > 0,
        lbl1(iw) = ic*ones(liw,1);
        rv = 1;
      end
    end
  end
  k = k + 1;
end

% Reorder cluster labels

[lbl1,ncl1] = ord_clus(lbl1);
