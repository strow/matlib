function r = nbresp(v, L, t)

% function r = nbresp(v, L, t)
%
% Norton-Beer response function, via cosine transform
% of the corresponding apodization function
%
% inputs
%   v - wavenumbers (equally spaced, ascending from 0)
%   L - max path length
%   t - weight selector (1-8, from nbapod.m)
%
% output
%   r - Norton-Beer response function of v

if nargin == 2
  t = 3;
elseif nargin == 1
  t = 3;
  L = 1;
end

L1 = L;
v = v(:);
n = length(v);

cospts = 2^nextpow2(n) + 1;
dv = v(2) - v(1);
vmax = dv * (cospts-1);
dd = 1/(2*vmax);
L1ind = round(L1/dd) + 1;
L1a = dd * (L1ind-1);
L1pts = 0:dd:L1a;

intf = zeros(cospts,1);
intf(1:L1ind) = nbapod(L1pts', L1, t);
spec = real(fft([intf; flipud(intf(2:cospts-1,1))]));
spec = spec(1:cospts) / max(spec(1:cospts));

r = spec(1:n);

