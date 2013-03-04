function c = kbresp(v, L, k)

% function c = kbresp(v, L, k)
%
% Kaiser-Bessel response function
%
% inputs
%   v - wavenumber; may be a vector
%   L - max path length
%   k - Kaiser-Bessel parameter (int > 0)
%
% output
%   c - Kaiser-Bessel response of v

if nargin == 2
  k = 6;
elseif nargin == 1
  k = 6;
  L = 1;
end

v = v(:);
n = length(v);

y = 2*pi*v*L;

t = 1 - (y/k).^2;

c = zeros(n,1);

for j = 1:n
  
  if t(j) > 0

    c(j) = sinh(k * sqrt(t(j))) / (sinh(k) * sqrt(t(j))) ;

  elseif t(j) < 0

    c(j) = sin(k * sqrt(-t(j))) / (sinh(k) * sqrt(-t(j))) ;

  else

    c(j) = 0;

  end
end

