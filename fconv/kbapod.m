function c = kbapod(d, L, k)

% function c = kbapod(d, L, k)
%
% Kaiser-Bessel apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%   k - Kaiser-Bessel parameter (int > 0)
%
% output
%   c - apodization of d

if nargin == 2
  k = 6;
elseif nargin == 1
  k = 6;
  L = 1;
end

[m,n] = size(d);

x = k * sqrt(1 - (d/L).^2) ;

% I(x)
r = ones(m,n); f = 1;
for j = 1:8;
  f = f * j;
  r = r + (x/2).^(2*j) / f^2 ;
  end

% I(k)
s = 1; f = 1;
for j = 1:8;
  f = f * j;
  s = s + (k/2).^(2*j) / f^2 ;
  end

c = (abs(d) <= L) .* r / s;

