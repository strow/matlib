%
% NAME
%   fconv1 - basic apodization 
% 
% SYNOPSIS
%   function [rad2, frq2] = fconv1(rad1, frq1, atype, opts)
%
% INPUTS
%    rad1   - channel radiances, in column order
%    frq1   - channel centers (vector of wavenumbers)
%    atype  - apodization type (a string)
%    opts   - optional interferometric and apodization param's 
%
% OUTPUTS
%    rad2   - convolved channel radiances
%    frq2   - channel center wavenumbers
%
% NOTES
%   fconv1 calls apod.m to return the desired apodization.
%   The input and output frequency grids are identical, and 
%   we assume dv = frq(2) - frq(1) evenly divides the grid 
%   points.  The convolution is performed as follows:
%
%         rad1 --> intf1 --> intf2 --> rad2 
%             ifft      apod       fft   
%
% AUTHOR
%   H. Motteler, 22 June 08

function [rad2, frq2] = fconv1(rad1, frq1, atype, opts)

% set defaults
dv = frq1(2) - frq1(1);
L = 1 / (2*dv);
aparg = 0;

% option to override defaults with opts fields
if nargin == 4
  optvar = fieldnames(opts);
  for i = 1 : length(optvar)
    vname = optvar{i};
    if exist(vname, 'var')
      % fprintf(1, 'fconv1: setting %s\n', vname)
      eval(sprintf('%s = opts.%s;', vname, vname));
    else
      fprintf(1, 'fconv1: WARNING unknown option variable %s\n', vname);
    end
  end
end

% guarantee frq1 is a column vector
frq1 = frq1(:);

% check that input args conform
[nrow,ncol] = size(rad1);
if nrow ~= length(frq1)
  error('rad1 rows must match frq1 length')
end

% get v1 and v2, and their indices
v1 = frq1(1);
v1ind = round(v1/dv);
v2 = frq1(nrow);
v2ind = round(v2/dv);

% find vmax and dd
vmaxind = 2^nextpow2(v2ind - 1);
vmax = (vmaxind - 1) * dv;
dd = 1 / (2*vmax);

% set the output frequencies
frq2 = frq1;

% initialize output array
rad2 = zeros(nrow, ncol);

% build the apodization vector
Lpts = (0:vmaxind-1)*dd;
apvec = apod(Lpts, L, atype, aparg)';

intf2 = zeros(vmaxind, 1);
rad1b = zeros(vmaxind, 1);

% loop on columns of rad1
for i = 1 : ncol

  % embed input data in a [0-vmax] interval
  rad1b(v1ind:v2ind) = rad1(:,i);

  % transform spectral data to an interferogram
  intf1 = real(ifft([rad1b; flipud(rad1b(2:vmaxind - 1,1))]));

  % apply the apodization
  intf2(1:vmaxind) = intf1(1:vmaxind) .* apvec;

  % transform the new interferogram back to radiances
  rad2a = real(fft([intf2; flipud(intf2(2:vmaxind-1,1))]));

  % save radiances from the input band
  rad2(:,i) = rad2a(v1ind:v2ind);

end

