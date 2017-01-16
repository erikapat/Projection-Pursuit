function [lbl1,ncl1] = ord_clus(lbl)

%
% [label_r,nclus_r] = ord_clus(label)
%
% Cluster identification from projections onto directions
% maximizing and minimizing the kurtosis coefficient of
% the data
%
% Housecleaning of cluster labels (to be used as
% a subroutine of clus_kur)
%
% Input:    label,    labels assigned to the observations
% Outputs:  label_r,  revised labels
%           nclus_r,  revised number of clusters
%

% DP/FJP  6/29/01

% Initialize parameters

n = length(lbl);
lbl1 = zeros(n,1);

% Order labels

[lbls,ix] = sort(lbl);

% Revise labels

r = 1;
s = 0;
t0 = 0;
while r <= n,
  t = lbls(r);
  if t ~= t0,
    s = s + 1;
    t0 = t;
  end
  lbl1(ix(r)) = s;
  r = r + 1;
end

% Number of different labels

ncl1 = s;
