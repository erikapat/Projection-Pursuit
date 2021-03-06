function [lblf,ncl] = clus_kur(x)

%
% [label,ncl] = clus_kur(x)
%
% Cluster identification from projections onto directions
% maximizing and minimizing the kurtosis coefficient of
% the data
%
% Master program
%
% Input:   Observations x, matrix with one observation in each row
% Output:  label, a numerical label assigning each observation
%                 to one cluster
%          ncl, number of identified clusters
%

% DP/FJP  6/29/01

% initial data

[n,p] = size(x);
lbl = ones(n,1);

% main cluster identification routine

lbl = clus_bas(x,1);

% study possible mergers of identified clusters

[lblf,ncl] = clus_grp(x,lbl);
