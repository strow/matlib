function b = cosapod(a, L)

% function b = cosapod(d, L)
%
% cosine apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%
% output
%   b - apodization of d
%
% apodization with cosine function from 0 to PI, normalized to [0,1]

if nargin == 1
  L = 1;
end

b = (abs(a) <= L) .* (1 + cos(pi*a/L)) / 2;

