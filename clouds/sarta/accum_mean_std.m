function [sumy,sumysqr,N] = accum_mean_std(sumx,sumxsqr,M,xn,ncount)

%% since I keep running out of memory when doing 25 x 2378 x 12150 runs for SARTA 100 layer, maximum random overlap,
%% this is a way of just keeping three matrices going, to allow me to do the std dev and mean calcs more efficiently
%% after which
%% mu = sumy./N; stddev = sqrt((sumysqr - 2*mu.*sumy + N.*mu.*mu)./(N-1));
%%
%% eg
%% [sumy,sumysqr,Nmatr] = accum_mean_std(0,0,0,p.rcalc,1);                             
%% [sumy,sumysqr,Nmatr] = accum_mean_std(sumy,sumysqr,Nmatr,p.rcalc+0.25*randn(size(p.rcalc)),2);
%% [sumy,sumysqr,Nmatr] = accum_mean_std(sumy,sumysqr,Nmatr,p.rcalc+0.25*randn(size(p.rcalc)),3);
%% [sumy,sumysqr,Nmatr] = accum_mean_std(sumy,sumysqr,Nmatr,p.rcalc+0.25*randn(size(p.rcalc)),4);
%% mu = sumy./Nmatr; sss = sqrt((sumysqr - 2*mu.*sumy + Nmatr.*mu.*mu)./(Nmatr-1));

%%
%% input
%%   x      = accumulated sum   (2378 x 12150) for past ncount-1 iterations
%%   xsqr   = accumulated sum^2 (2378 x 12150) for past ncount-1 iterations
%%   xn     = current (2378 x 12150) we want to include at this ncount th iteration
%%   ncount = current iteration
%%   M      = of the 2378 x 12150 channels, sometimes they are NaN or Inf and so should not be
%%            included in the averaging; M is a count of how many are actauuly included
%%
%% output
%% y        = x + xn
%% ysrq     = xsqr + xn^2
%% N        = M + 1 (in a matrix sense)

if ncount == 1
  sumy    = zeros(size(xn));
  sumysqr = zeros(size(xn));
  N       = zeros(size(xn));

  boo          = find(isfinite(xn));
  sumy(boo)    = xn(boo);
  sumysqr(boo) = xn(boo).*xn(boo);
  N(boo)       = 1;

else
  sumy    = sumx;
  sumysqr = sumxsqr;
  N       = M;

  boo          = find(isfinite(xn));
  sumy(boo)    = xn(boo) + sumx(boo);
  sumysqr(boo) = xn(boo).*xn(boo) + sumxsqr(boo);
  N(boo)       = M(boo) + 1;

  %[xn(1,1) sumy(1,1) sumysqr(1,1) N(1,1)]
end
