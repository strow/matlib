function [rad2, wnum2] = fconv(rad1, wnum1, ifp, atype, aparg)

% function [rad2, wnum2] = fconv(rad1, wnum1, ifp, atype, aparg)
%
% Fourier convolution of channel radiances
%
% Inputs
%    rad1   - channel radiances, in column order
%    wnum1  - channel centers (vector of wavenumbers)
%    ifp    - interferometric parameter file (default "ifp.m")
%    atype  - optional apodization spec (default Kaiser-Bessel)
%    aparg  - optional apodization parameter 
%
% Outputs
%    rad2   - convolved channel radiances
%    wnum2  - channel center wavenumbers
%
% The convolution is performed as follows:
%
%   rad1 --> intf1 --> intf2 --> rad2 
%       ifft      apod       fft   
%
% fconv() is intended for the case where the input and output
% channel sets are either the same, or where the input channel 
% set is a subset of the output set.
% 
% The channel set (passed in as wnum1) does not have to be exactly
% the same as the [v1,v2] band of the parameter file; in fact it
% can be any monotonically increasing set of multiples of dvc less
% than vmax.
%
% See apod.m for available apodizations.  If no apodization type 
% or parameter is specified, Kaiser-Bessel #6 is used.  Most of
% the apodizations do not need the parameter aparg.

% H. Motteler, 8/24/98

% warning off by LLS, tired of warnings filling screen
warning off

% get interferometric parameters
%
calcifp(ifp)
global L1 L2 L1ind L2ind dvc dd

% set input defaults
%
if nargin == 4
  aparg = -1;
elseif nargin == 3
  atype = 'kb';
  aparg = 6;
elseif nargin == 2
  atype = 'kb';
  aparg = 6;
  ifp   = 'ifp.m';
end  

% get channel indices relative to 0:dvc:vcmax
%
vcmaxind = 2^nextpow2(L2ind-1) + 1;
vcmax = (vcmaxind - 1) * dvc;
wnum0 = 0:dvc:vcmax;

chan_ind = interp1(wnum0, 1:vcmaxind, wnum1, '*nearest');
% chan_ind = unique(chan_ind(~isnan(chan_ind)));
wnum2 = wnum0(chan_ind);
nchans = length(wnum2);

% initialize arrays
%
rad0  = zeros(vcmaxind, 1);
intf2 = zeros(vcmaxind, 1);

[m, n] = size(rad1);
rad2 = zeros(nchans, n);

% build the apodization array
%
L2pts = 0:dd:(L2ind-1)*dd;
apvec  = apod(L2pts, L2, atype, aparg)';


% convolve each column of rad1:
%
% rad1 --> intf1 --> intf2 --> rad2 
%     ifft      apod       fft   
% 
for i = 1:n

  rad0(chan_ind) = rad1(:,i);

  intf1 = real(ifft([rad0; flipud(rad0(2:vcmaxind-1,1))]));

  intf2(1:L2ind) = intf1(1:L2ind) .* apvec;

  rad3 = real(fft([intf2; flipud(intf2(2:vcmaxind-1,1))]));

  rad2(:,i) = rad3(chan_ind);

end

