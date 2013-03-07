%
% NAME
%   bt2rad - translate brightness temperature to radiance
%
% SYNOPSIS
%   rad = bt2rad(fr, bt);
%
% INPUTS
%   fr  - n-vector of wavenumbers, cm-1
%   bt  - n x k array of brightness temps, K
%
% OUTPUT
%   rad - n x k array of radiances, mW/m2 per strad
%
% DISCUSSION
%   radiance units are milliwatts per square meter per steradian.
%
%   fr can be a row or column vector, and bt a vector or array.
%   When the length of fr equals the number of rows of bt, it is
%   applied along columns.  When fr is scalar it is applied to all
%   elements of bt.  If bt is a row vector with the same length 
%   as fr, then for this one case fr is applied along the row.
%
% AUTHOR
%   H. Motteler
%

function rad = bt2rad(fr, bt);

% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1

% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;
c2 = (h*c/k) * 100;

% the scalar calculation, for reference
% rad = c1 * fr^3 / (exp(c2*fr/bt) - 1)

% set up the data
fr = fr(:);               % make fr a column vector
k = length(fr);           % fr vector length
dd = size(bt);            % bt dimension list
[m,n] = size(bt);         % bt size as a 2d array
bt = reshape(bt, m, n);   % make bt a 2d array

% special case: bt is a row vector and length(bt) = length(fr)
if m == 1 && n == k
  bt = bt'; n = 1; m = k;  
end

% special case: fr is a scalar and bt has more than one row
if k == 1 && m > 1
  fr = fr * ones(m,1); k = m;
end

% check that fr conforms with bt
if m ~= k
  error('the length of fr must equal the number of rows of bt')
end

% do the vectorized calculation
rad  = c1 * (fr.^3 * ones(1,n)) ./ (exp((c2 * fr * ones(1,n)) ./ bt) - 1);

% restore bt original shape
rad = reshape(rad, dd);

