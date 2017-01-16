function [y,m,R] = normaliz(x)

%
% [y,m,R] = normaliz(x)
%
% Se generan datos con media cero y matriz de covarianza la identidad
% Los datos de partida vendran dados como filas de una matriz
%

[n,p] = size(x);

if n < p,
  disp('Mas dimensiones que datos');
  break
end

m = mean(x);
y = x - ones(n,1)*m;
S = cov(y);
R = chol(S);
y = ((R')\(y'))';
