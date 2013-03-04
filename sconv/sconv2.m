
function [rout, fout] = sconv2(rin, fin, clist, sfile, coff);

% NAME
%
%   sconv2 -- do convolutions directly from tabulated SRFs
%
% SYNOPSIS
%
%   function [rout, fout] = sconv2(rin, fin, clist, sfile, coff);
%
% INPUTS
%
%   rin    - m x n input radiance data
%   fin    - m-vector of input frequencies
%   clist  - k-vector of output channel numbers
%   sfile  - filename of HDF-format SRF data
%   coff   - optional channel center offsets
% 
% OUTPUTS
%
%   rout   - k x n convolved ouput data
%   fout   - k-vector of output frequencies
% 
% DESCRIPTION
% 
%   sconv2 calculates an SRF convolution by interpolating and
%   calculating inner products from  the tabulated SRFs directly.
%   It performs the same function as sconv1, without needing to 
%   first read a sparse convolution matrix; it uses less memory
%   but is a little slower for more than a few hundred channels.
%
% H. Motteler, 15 Jan 02
%

% guarantee input vectors are columns
fin = fin(:);
clist = clist(:);

% check rin dimensions
[m,n] = size(rin);
if m ~= length(fin)
  error('number of rows of rin do not match length of fin')
end

% set a default HDF SRF file
if nargin < 4
  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
end

% default to zero channel offsets
if nargin < 5
  coff = zeros(length(clist), 1);
else
  coff = coff(:);
  if length(coff) ~= length(clist)
    error('coff and clist must be the vectors of the same length');
  end
end

% read the SRF data
[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   freq   = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end
clear alist fattr

% use the subset of channels specified in clist
cind = interp1(chanid, 1:length(chanid), clist, 'nearest');
chanid = chanid(cind);
freq = freq(cind) + coff;
srfval = srfval(cind,:);
width = width(cind);
[nchan, nspts] = size(srfval);

% srfreq is an nchan x nspts array of frequency points for srfval
srfreq = (width * fwgrid) + freq * ones(1, nspts);

% set the output frequency list
fout = freq;

% check that fin is uniformly spaced and increasing
dfin = diff(fin);
dmax = max(dfin);
dmin = min(dfin);
if abs(dmax - dmin) > 10e-9
  error('fin must be uniformly spaced')
end

% get the step size of fin
dv = dfin(1);
if dv <= 0
  error('fin must be in increasing order')
end

% get the span of fin
nfin = length(fin);
f1 = fin(1);
f2 = fin(nfin);

% initialize output array
[m,n] = size(rin);
rout = zeros(nchan, n);

% loop on SRFs in cind
for j = 1 : nchan

  % get the frequency span of the current SRF
  v1 = srfreq(j, 1);
  v2 = srfreq(j, nspts);

  % if the SRF is outside fin, skip this channel
  if v2 <= f1 | f2 <= v1
    fprintf(1, 'sconv2(): WARNING -- SRF %d outside of input range\n', ...
	    chanid(j));
    continue
  end

  % if the SRF overlaps fin, lop it off to fit; 
  % giv a warning message if we lop more than dv
  if v1 < f1
    if v1 < f1 - dv
      fprintf(1, 'sconv2(): WARNING -- truncating LHS of SRF %d\n', ...
	      chanid(j));
    end
    v1 = f1;
  end
  if f2 < v2
    if f2 + dv < v2
      fprintf(1, 'sconv2(): WARNING -- truncating RHS of SRF %d\n', ...
	      chanid(j));
    end
    v2 = f2;
  end

  % find the indices of the current SRF in fin 
  v1ind = ceil((v1-f1)/dv) + 1;
  v2ind = floor((v2-f1)/dv) + 1;

  % interpolate the SRF to a subinterval of the fin grid
% s1 = interp1(srfreq(j,:), srfval(j,:), fin(v1ind:v2ind), 'spline');
  s1 = interp1(srfreq(j,:), srfval(j,:), fin(v1ind:v2ind), 'linear');

  % normalize the SRF
  s1 = s1 ./ sum(s1);

  % apply the SRF
  rout(j, :) = s1' * rin(v1ind:v2ind, :);
end

