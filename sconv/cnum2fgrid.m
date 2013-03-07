
function fgrid = cnum2fgrid(clist, sfile);

% NAME
%
%   cnum2fgrid  -- get a 0.0025 1/cm grid from a channel list
%
% SYNOPSIS
%
%   function fgrid = cnum2fgrid(clist, sfile);
%
% INPUTS
%
%   clist  - channel number list 
%   sfile  - HDF SRF tabulation data 
% 
% OUTPUTS
%
%   fgrid   - vector of output frequencies
% 
% DESCRIPTION
% 
%   cnum2fgrid takes a list of channel numbers from an SRF 
%   tabulation and returns a 0.0025 1/cm grid bounded above
%   and below by the listed SRFs, including the SRF wings.
%
% H. Motteler, 15 Jan 02
%

% the default channel list is all 2378 AIRS channels
if nargin < 1
  clist = 1 : 2378;
end

% make clist a column vector
clist = clist(:);

% the default SRF file is the V1.0 tabulation
if nargin < 2
  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
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
freq = freq(cind);
srfval = srfval(cind,:);
width = width(cind);
[nchan, nspts] = size(srfval);

% srfreq is an nchan x nspts array of frequency points for srfval
srfreq = (width * fwgrid) + freq * ones(1, nspts);

% get spanning frequencies for this channel set
vmin = min(min(srfreq));
vmax = max(max(srfreq));

% calculate the spanning band on a .0025 1/cm grid
dv = 0.0025;
v1 = ceil(vmin/dv)*dv;
v2 = floor(vmax/dv)*dv;
fgrid = v1 : dv : v2;

