function [rch, wch] = fconvkc(rkc, ifp, atype, aparg)

% function [rch, wch] = fconvkc(rkc, ifp, atype, aparg)
% 
% Fourier convolution of kcarta monochromatic radiances
%
% inputs
%   rkc   - kcarta monochromatic radiances, in column order
%   ifp   - interferometric parameter file (default "params.m")
%   atype - optional apodization type (default Kaiser-Bessel)
%   aparg - optional apodization parameter
%
% outputs
%   rch   - convolved channel radiances
%   wch   - channel center wavenumbers
%
% The convolution is performed as follows:
%
%   rkc --> rVL --> intf1 --> intf2 --> rch
%     interp   ifft      apod       fft   
%
% Convolution parameters, including the monochromatic interval
% [v1,v2], are taken from the interferometric parameter file. 
% It is assumed that rkc will span the entire interval.
% 
% See apod.m for available apodizations.  If no apodization type 
% or parameter is specified, Kaiser-Bessel #6 is used.  Most of
% the apodizations do not need the parameter aparg.
%
% This version does a linear interpolation from kcarta to vlaser 
% point spacings.

% set input defaults
%
if nargin == 3
  aparg = -1;
elseif nargin == 2
  atype = 'kb';
  aparg = 6;
elseif nargin == 1
  atype = 'kb';
  aparg = 6;
  ifp   = 'params';
end  

% compute interferometric parameters
%
calcifp(ifp);
global v1 v2 vmax rolloff
global L1 L2 Lcut Lmax
global dvk dv dvc ddk dd
global v1ind v2ind vmaxind vLpts
global L1ind L2ind Lcutind Lmaxind
global v1cind v2cind nchans
global v1kind v2kind kcpts

% check supplied radiance data
%
[m,n] = size(rkc);
if m ~= kcpts
  error('input rows do not match interferometric parameters')
end

% frequency domain filter at band edges
%
rpts = round(rolloff / dv);  % number of rolloff points
f0 = (1+cos(pi+(0:rpts-1)*pi/(rpts)))/2;
filt = [f0, ones(1,v2ind-v1ind+1-2*rpts), fliplr(f0)];

% initialize arrays
%
wkc = dvk*((v1kind:v2kind)-1);

r1vL  = zeros(vLpts, 1); 
r2vL  = zeros(vmaxind, 1);
wvL =  dv*((v1ind:v2ind)-1);

% check edges of interpolation interval
if wvL(1) < wkc(1)
  wvL(1) = wkc(1);
end
if wkc(kcpts) < wvL(vLpts)
  wvL(vLpts) = wkc(kcpts);
end

rch = zeros(nchans, n);
wch = dvc*((v1cind:v2cind)-1);

% build the apodization array
%
L1pts = 0:dd:(L1ind-1)*dd;
apvec  = apod(L1pts, L1, atype, aparg)';
intf2 = zeros(L2ind, 1);	

% convolve each column of rkc:
%
% rkc --> rVL --> intf1 --> intf2 --> rch
%   interp   ifft      apod       fft   
% 
for i = 1:n

  r1vL = interp1(wkc, rkc(:,i), wvL, '*linear');

  r2vL(v1ind:v2ind) = r1vL .* filt;

  intf1 = real(ifft([r2vL; flipud(r2vL(2:vmaxind-1,1))]));

  intf2(1:L1ind) = intf1(1:L1ind) .* apvec;

  rtmp = real(fft([intf2; flipud(intf2(2:L2ind-1,1))]));

  rch(:,i) = rtmp(v1cind:v2cind);

end

