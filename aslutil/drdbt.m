function [deriv]=drdbt(freq,bt);
% function [deriv]=drdbt(freq,bt);
%    dR/dB(T), used by shiftdbt.m
c1=1.1911E-8;
c2=1.4387863;
deriv=c1.*c2.*freq.^4.*exp(c2.*freq./bt)./(bt.^2.*(exp(c2.*freq./bt)-1.0).^2);
