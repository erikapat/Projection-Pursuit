%
% Script to generate the simulation results in the paper
% "Outlier Identification using Projections," by D. Peña
% and F.J. Prieto
%

% DP/FJP 23/5/00

p = [ 4 8 15 30 ];
k = [ 2 4 8 ];
dst0 = [ 7 10 14 ; 6 9 13 ; 5 8 12 ; 4 7 11 ];
%nrep = 100;
nrep = 5;
rslt = [];
outp = 'res0.txt';

for md = 0:3,
  mode = -md;
  eval(['diary ' outp]);
  if md == 0,
    disp('Normal observations w/ outliers - kurtosis');
  elseif md == 1,
    disp('Normal observations - kurtosis');
  elseif md == 2,
    disp('Uniform observations - kurtosis');
  else
    disp('Student-t observations - kurtosis');
  end
  disp(' Num   Dim  Clus   Dist   Ediv   Emix   Eboth      Nerr    Rmax   Rmean');
  diary off
  dst = 2*dst0;
  rslt = [];
  for i0=1:length(p),
    pp = p(i0);
    for i1=1:length(k),
      kk = k(i1);
      dd = dst(i0,i1);
      nn = 20*pp;
      v = sim_clus(nn,pp,kk,dd,nrep,mode,0.1);
      w = [ nn pp kk dd v ];
      eval(['diary ' outp]);
      aa = sprintf(' %3.0f   %3.0f  %3.0f ',nn,pp,kk);
      bb = sprintf(' %6.2f  %5.2f  %5.2f ',dd,v(1),v(2));
      cc = sprintf('  %5.2f  %8.2f %7.2f ',v(3),v(4),v(5));
      dd = sprintf('%7.2f',v(6));
      disp([ aa bb cc dd ]);
      diary off
      rslt = [ rslt ; w ];
    end
  end
end

for md = 0:3,
  mode = -md;
  eval(['diary ' outp]);
  if md == 0,
    disp('Normal observations w/ outliers - Jones/Sibson');
  elseif md == 1,
    disp('Normal observations - Jones/Sibson');
  elseif md == 2,
    disp('Uniform observations - Jones/Sibson');
  else
    disp('Student-t observations - Jones/Sibson');
  end
  disp('Num   Dim   Clus   Dist   Ediv   Emix   Eboth   Nerr   Rmax   Rmean');
  diary off
  dst = 2*dst0;
  rslt = [];
  for i0=1:length(p),
    pp = p(i0);
    for i1=1:length(k),
      kk = k(i1);
      dd = dst(i0,i1);
      nn = 20*pp;
      v = sim_cljs(nn,pp,kk,dd,nrep,mode,0.1);
      w = [ nn pp kk dd v ];
      eval(['diary ' outp]);
      aa = sprintf(' %3.0f   %3.0f  %3.0f ',nn,pp,kk);
      bb = sprintf(' %6.2f  %5.2f  %5.2f ',dd,v(1),v(2));
      cc = sprintf('  %5.2f  %8.2f %7.2f ',v(3),v(4),v(5));
      dd = sprintf('%7.2f',v(6));
      disp([ aa bb cc dd ]);
      diary off
      rslt = [ rslt ; w ];
    end
  end
end
