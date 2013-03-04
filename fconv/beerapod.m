function b = beerapod(d, L)

% function b = beerapod(d, L)
%
% Beer apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%
% output
%   b - apodization of d
% 
% NOTE: from the Wisconson group, via Larrabee

b = (abs(d) <= L) .* (1 - (d/L).^2).^2;

