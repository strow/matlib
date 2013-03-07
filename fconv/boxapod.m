function y = boxapod(d, L)

% function y = boxapod(d, L)
%
% boxcar apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%
% output
%   y - apodization of d


if nargin == 1
  L = 1;
end

%y = (abs(d) <= L);

y = double( (abs(d) <= L) );
