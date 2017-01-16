function [v,lam] = max_js(x,d0)

%
% [v,lambda] = max_js(x,d)
%
% Cluster identification from projections onto directions
% maximizing the Jones and Sibson criterion
%
% Computation of one projection direction
%
% Inputs:   observations, x, matrix with one observation in each row
%           d, initial estimate of the direction (optional)
% Output:   v, direction corresponding to a local maximizer of the
%              criterion
%           lambda, multiplier for the normalization constraint
%

% DP/FJP  6/29/01

% Initial values

if nargin < 3,
  show = 0;
end

if nargin < 2,
  d0 = [];
end

maxit = 100;
tol = 1.0e-5;
tol1 = 1.0e-6;
beta = 1.0e-4;
rho0 = 0.1;

[n,p] = size(x);
mdi = 1 + floor(n/2);

mm = mean(x);
S = cov(x);
xx = x;
if norm(mm) > 1000*eps,
  xx = x - ones(n,1)*mm;
end
R = eye(p);
if norm(R - S) > 1.0e4*eps,
  R = chol(S);
  xx = xx*inv(R);
end

if length(d0) == 0,
  uu = zeros(n,p);
  for j = 1:n,
    uu(j,:) = xx(j,:)/(eps + norm(xx(j,:)));
  end
  Su = cov(uu);
  [V,D] = eig(Su);
  [v,ik] = max(diag(D));
  a = V(:,ik);
else
  a = R*d0;
  a = a/norm(a);
end

z = xx*a;
sk3 = sum(z.^3);
sk4 = sum(z.^4) - 3;
sk = sk3^2 + 0.25*sk4^2;

al = 0;
it = 0;
diff = 1;
rho = rho0;
clkr = 0;

% Iteration of the optimization procedure

if show,
  disp(' It.      F.obj.       | g |           c         alfa       rho');
end

while 1,

  nmd2 = a'*a;
  c = nmd2 - 1;

  gaux = (6*sk3 + 2*sk4*z).*(z.^2);
  g = xx'*gaux;
  lam = 0.5*g'*a;
  f = sk - lam*c - 0.5*rho*c^2;
  gl = g - 2*lam*a;

  if show,
    aa = sprintf('%3.0f  %12.5f %13.4e',it,f,norm(gl));
    bb = sprintf(' %13.4e %8.3f %11.2e',abs(a'*a-1),al,rho);
    disp([ aa bb ]);
  end

  A = 2*a';
  [Q,W] = qr(a);
  Z = Q(:,(2:p));

  crit = norm(gl) + abs(c);
  if (crit <= tol)|(it >= maxit),
    break
  end

  aux1 = (z.^2)'*xx;
  aux0 = 18*aux1'*aux1;
  aux2 = (z.^3)'*xx;
  aux0 = aux0 + 8*aux2'*aux2;
  aux3 = (((12*sk3 + 6*sk4*z).*z)*ones(1,p)).*xx;
  H = aux0 + xx'*aux3;
  Hl = H - 2*lam*eye(p,p);
  Hr = Z'*Hl*Z;
  [V,E] = eig(Hr);
  Es = min(-abs(E),-1.0e-4);
  Hs = V*Es*V';

  py = - c/(A*A');
  rhs = Z'*(g + H*A'*py);
  pz = - Hs\rhs;
  pp = Z*pz + py*A';

  dlam = (2*a)\(gl + H*pp);

  al = 1;
  f0d = gl'*pp - 2*rho*c*a'*pp - dlam*c;
  crit1 = beta*norm(pp)^2;

  if f0d < crit1,
    rho1 = 2*(crit1 - f0d)/(eps + c^2);
    rho = max([2*rho1 1.5*rho rho0]);
    f = sk - lam*c - 0.5*rho*c^2;
    f0d = gl'*pp - 2*rho*c*a'*pp - dlam*c;
    clkr = 0;
  elseif (f0d > 1000*crit1)&(rho > rho0),
    rho1 = 2*(crit1 - gl'*pp + dlam*c)/(eps + c^2);
    if (clkr == 4)&(rho > 2*rho1),
      rho = 0.5*rho;
      f = sk - lam*c - 0.5*rho*c^2;
      f0d = gl'*pp - 2*rho*c*a'*pp - dlam*c;
      clkr = 0;
    else
      clkr = clkr + 1;
    end
  end
  if (abs(f0d) < tol1),
    break
  end

  itbl = 0;
  while itbl < 20,
    aa = a + al*pp;
    lama = lam + al*dlam;
    zz = xx*aa;
    cc = aa'*aa - 1;
    sk3 = sum(zz.^3);
    sk4 = sum(zz.^4) - 3;
    skk = sk3^2 + 0.25*sk4^2;
    ff = skk - lama*cc - 0.5*rho*cc^2;
    if ff > f + 0.01*al*f0d,
      break
    end
    al = al/2;
    itbl = itbl + 1;
  end
  if itbl >= 20,
    disp('Failure in the line search');
    break
  end

  a = aa;
  lam = lama;

  z = xx*a;
  sk3 = sum(z.^3);
  sk4 = sum(z.^4) - 3;
  sk = sk3^2 + 0.25*sk4^2;

  it = it + 1;

end

a = inv(R)*a;
v = a/norm(a);