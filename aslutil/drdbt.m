function [deriv]=drdbt(freq,bt);
% function [deriv]=drdbt(freq,bt);
%    dR/dB(T), used by shiftdbt.m

% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1

% Compute radiation constants c1 and c2
% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;
c2 = (h*c/k) * 100;

bt = bt(:);
freq = freq(:);
deriv=c1.*c2.*freq.^4.*exp(c2.*freq./bt)./(bt.^2.*(exp(c2.*freq./bt)-1.0).^2);
