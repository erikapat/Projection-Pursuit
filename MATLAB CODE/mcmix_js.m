function [T,Vv] = mcmix_js(x0)

%
% [T,V] = mcmix_js(x)
%
% Cluster identification from projections onto directions
% maximizing the Jones and Sibson criterion
%
% Computation of projection directions maximizing the criterion
%
% Inputs:   observations, x, matrix with one observation in each row
% Outputs:  T, projections of each observation onto each projection
%              direction
%           V, projection directions maximizing the criterion
%

% DP/FJP  6/29/01

maxit = 100;

% Initialization

tol = 1.0e-5;
tol1 = 1.0e-6;
beta = 1.0e-4;
rho0 = 0.1;

[n,p0] = size(x0);
mdi = 1 + floor(n/2);

T = zeros(n,2*p0);
ti = 1;

x = x0;
p = p0;
pin = p - 1;

mm = mean(x);
S = cov(x);
xx = x - ones(n,1)*mm;
R = chol(S);
xx = xx*inv(R);

for i = 1:pin,

  a = max_js(xx);
  [Q,W] = qr(a);

% Outlyingness of projections

  Tt = xx*Q;

  T(:,ti) = Tt(:,1);
  eval(['Q' num2str(ti) ' = Q;']);
  ti = ti + 1;

% Reduction of dimension

  xx = Tt(:,(2:p));
  p = p - 1;

end

% Generation of the last projection

T(:,ti) = xx;
ti = ti + 1;

% Regenerating information on the projections

ti = p0 - 1;

eval(['U = Q' num2str(ti) ';']);

ti = ti - 1;

while ti > 0,

  tk = p0 - ti + 1;

  eval(['Z = Q' num2str(ti) '(:,(2:tk));']);
  eval(['a = Q' num2str(ti) '(:,1);']);
  U1 = Z*U;
  U = [ a U1 ];

  ti = ti - 1;

end

U = inv(R)*U;
uaux = diag(U'*U);
Vv = U*diag(1../sqrt(uaux));
