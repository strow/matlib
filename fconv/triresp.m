function b = triresp(v, L)

% function b = triresp(v, L)
%
% triangle response function
%
% inputs
%   v - wavenumber; may be a vector
%   L - max path length
%
% output
%   b - response function of v


if nargin == 1
  L = 1;
end

q = (v == 0);

v1 = v + q;

b = ~q .* (sin(pi*L*v1).^2 ./ (pi*L*v1).^2) + q;

