function a = nbapod(d, L, t)

% function a = nbapod(d, L, t)
%
% Norton-Beer apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%   t - weight selector (from table, below)
%
% output
%   a - apodization of d

if nargin == 2
  t = 3;
elseif nargin == 1
  t = 3;
  L = 1;
end

% tabulated Norton-Beer coefficients
%
c = zeros(8,5);
c(1,:) = [1         0         0         0        0       ]; % sinc
c(2,:) = [0.23977   0.45806   0.22498   0.07719  0       ]; % Connes
c(3,:) = [0.548    -0.0833    0.5353    0        0       ]; % F1, weak
c(4,:) = [0.384093 -0.087577  0.703484  0        0       ]; % F1, weak-2
c(5,:) = [0.26     -0.154838  0.894838  0        0       ]; % F2, medium
c(6,:) = [0.152442 -0.136176  0.983734  0        0       ]; % F2, medium-2
c(7,:) = [0.09      0         0.5875    0        0.32250 ]; % F3, strong
c(8,:) = [0.045335  0         0.554883  0        0.399782]; % F3, strong-2


a = zeros(size(d));

r = (d / L).^2 ;

for i = 1:5

  a = a + c(t,i) * (1-r).^(i-1) ;

end

