function Vv = kur_nwa(x,mode)

%
%    V = kur_nwa(x,mode)
%
% Cluster identification from projections onto directions
% maximizing and minimizing the kurtosis coefficient of
% the data (to be used as a subroutine of clus_kur)
%
% Computation of directions that maximize and minimize
% the kurtosis coefficient of the projections
% Data is assumed to be standardized in advance
%
% Inputs:  x, observations (by rows)
%          mode, controls the computation of only maximization
%                or maximization and minimization directions
% Outputs: V, directions maximizing/minimizing kurtosis coef.
%             maximization directions first
%

% DP/FJP 23/5/00

if nargin < 2,
  mode = 0;
end

% Initializations

%% Parameters (tolerances)

maxit = 100;

tol = 1.0e-5;
tol1 = 1.0e-6;
beta = 1.0e-4;
rho0 = 0.1;

[n,p0] = size(x);

%% Initialization of vectors

Vv = [];
ti = 1;

% Computing directions

%% Choice of minimization/maximization

cff = 1;
if mode < 0,
  cffmax = 2;
else
  cffmax = 3;
end

%% Main loop to compute 2p directions

while (cff < cffmax),

  xx = x;
  p = p0;
  pin = p - 1;
  M = eye(p0);

  for i = 1:pin,
    if cff == 1,
      a = max_kur(xx);
    else
      a = min_kur(xx);
    end
    la = length(a);
    za = zeros(la,1); za(1) = 1;
    w = a - za; nw = w'*a;
    if abs(nw) > eps,
      Q = eye(la) - w*w'/nw;
    else
      Q = eye(la);
    end

%% Compute projected values

    Vv = [ Vv (M*a) ];

    Qp = Q(:,2:p);
    M = M*Qp;
    ti = ti + 1;

%% Reduce dimension

    Tt = xx*Q;
    xx = Tt(:,(2:p));

    p = p - 1;

  end

%% Compute last projection

  Vv = [ Vv M ];
  ti = ti + 1;

%% Proceed to minimization

  cff = cff + 1;

end
