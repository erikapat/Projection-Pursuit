function mrk = clus_bas1(x0,show,cutoff)

%
% lbl = clus_bas(x,show,cutoff)
%
% Cluster identification from projections onto directions
% maximizing and minimizing the kurtosis coefficient of
% the data
%
% Computation of projections and analysis of the gaps (to be used as
% a subroutine of clus_kur)
%
% Inputs:   observations, x, matrix with one observation in each row
%           show, if = 1 output is generated for each iteration
%           cutoff, cutoff level for the identification of clusters
%                   on the projections
% Output:   lbl, labels assigning each observation to one cluster
%

% DP/FJP  6/29/01

if nargin < 3,
  cutoff = 0.1;
end
if nargin < 2,
  show = 0;
end

% Initial data

maxit = 100;
tol = 1.0e-5;
[n,p] = size(x0);

% Check for sample size

mrk = ones(n,1);
obs_elim = p - 2;  % Do not consider extreme observations
if n <= (2*obs_elim),
  return
end

if show,
  lbl = ['Cluster identification procedure'];
  disp(lbl);
  disp('Number of observations / dimension');
  disp([ n p ]);
end

% criteria to check gaps in the projections

alphac = cutoff/(p^3.322);
delta = 1 - alphac^(1/n);

% Computation of projection directions

%% Standardize data

mm = mean(x0);
S = cov(x0);
x = x0 - ones(n,1)*mm;

Rr = chol(S);
x = ((Rr')\(x'))';

%% Invoke the routine to compute projection directions

V = kur_nwa(x);
[nt,mt] = size(V);

ref = 1;
k = 1;

% Analysis of the gaps in the projections
k=4;
while (k <= mt),

%% Compute projections

  t = x*V(:,k);

  [ts,ix] = sort(t);

%% inverse transformation to uniforms

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

%% Check significance levels for gaps and assign cluster labels

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
