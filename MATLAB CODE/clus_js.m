function [lblf,ncl] = clus_js(x0)

%
% [lbl,nclus] = clus_js(x)
%
% Cluster identification from projections onto directions
% maximizing the Jones and Sibson criterion
%
% Master program
%
% Input:    observations, x, matrix with one observation in each row
% Outputs:  lbl, labels assigning each observation to one cluster
%           nclus, number of identified clusters
%

% DP/FJP  6/29/01

if nargin < 2,
  show = 0;            % show graphical output in each iteration
end

% Initialization

maxit = 100;
tol = 1.0e-5;
[n,p] = size(x0);

% Remove cases with too few observations

mrk = ones(n,1);
obs_elim = p - 2;
if n <= (2*obs_elim),
  return
end

if show,
  lbl = ['Multivariate cluster detection'];
  disp(lbl);
  disp('Number of observations / dimension');
  disp([ n p ]);
end

% criteria for gap thresholds

cutoff = 0.1;
alphac = cutoff/(p^3.322);
delta = -log(alphac)/(n+1);

% Computing projection directions

%% Standardize data

mm = mean(x0);
S = cov(x0);
x = x0 - ones(n,1)*mm;

Rr = chol(S);
x = ((Rr')\(x'))';

[T,V] = mcmix_js(x);
[nt,mt] = size(V);

% Analysis of projections

ref = 1;
k = 1;

while (k <= mt),

% Projections and gaps

  t = x*V(:,k);
  [ts,ix] = sort(t);

%% Transformation to uniforms (assuming normality)

  rt = 0.5*(1 + erf(ts/sqrt(2)));
  rdif = diff(rt);
  idf = find(rdif > delta);

  if show,
    [vdif,idif] = max(rdif);
    disp(idif);
    disp([vdif delta]);
    figure(1); plot(t,'*');
    figure(2); plot(1:(n-1),rdif,'r-'); hold on;
    plot(1:(n-1),delta*ones(n-1,1),'b-');
    hold off
    pause
  end

  idf = [0 ; idf ; n ];
  mm1 = length(idf);

  if mm1 > 2,
    rr = 1;
    while (rr < (mm1-1)),
      wref = -1;
      lwl = idf(rr) + 1; upl = idf(rr+1);
      auxix = ix(lwl:upl);
      auxw = mrk(auxix);
      [auxw1,iy] = sort(auxw);
      for ss = 1:(upl-lwl+1),
        tt = auxix(iy(ss));
        vref = mrk(tt);
        if vref ~= wref,
          ref = ref + 1;
          wref = vref;
        end
        mrk(tt) = ref;
      end
      rr = rr + 1;
    end
  end

  k = k + 1;

end

[lblf,ncl] = clus_grp(x,mrk);
