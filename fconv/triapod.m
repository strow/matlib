function b = triapod(a, L)

% function b = triapod(a, L)
%
% Kaiser-Bessel apodization
%
% inputs
%   a - distance; may be a vector
%   L - max path length
%
% output
%   b - apodization of a


if nargin == 1
  L = 1;
end

b = (abs(a) <= L) .* (L - a) / L;

