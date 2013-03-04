function y = hamresp(v, L)

% function y = hamresp(v, L)
%
% Hamming response function
%
% inputs
%   v - wavenumber; may be a vector
%   L - max path length
%
% output
%   y - hamming response at v


if nargin == 1
  L = 1;
end

a = 0.23;
b = 1 - 2*a;
x = 2*pi*L*v;

q =  x == pi;

xsq = x.*x + q;           % x^2 everywhere except where x == pi

y = b * sinxx(x) .* (1 + 2*a*xsq ./ (b*(pi*pi-xsq)));

y = y .* (1-q) + a * q;   % a where x == pi

y = y / 0.54 ;		  % normalize

