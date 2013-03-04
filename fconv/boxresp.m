function r = boxresp(v, L)

% function r = boxresp(v, L)
%
% boxcar response function
%
% inputs
%   v - wavenumbers
%   L - max path length
%
% output
%   r - boxcar response at v


if nargin == 1
  L = 1;
end

r = sinxx(2*pi*L*v);

