
function [v1, v2] = cnum2vspan(clist, sfile);

% NAME
%
%   cnum2vspan -- get spanning frequency bounds from a channel list
% 
% SYNOPSIS
%
%   function [v1, v2] = cnum2vspan(clist, sfile);
%
% INPUTS
%
%    clist  - channel number list 
%    sfile  - HDF SRF tabulation data 
% 
% OUTPUTS
%
%    v1     - lower frequency bound on clist SRF set
%    v2     - upper frequency bound on clist SRF set
% 
% DESCRIPTION
% 
%    cnum2vspan takes a list of channel numbers from an SRF 
%    tabulation and returns lower and upper frequency bounds,
%    including the SRF wings, for the channels specified.
%
% H. Motteler, 15 Jan 02
%

% the default channel list is all 2378 AIRS channels
if nargin < 1
  clist = 1 : 2378;
end

% make clist a column vector
clist = clist(:);

% the default SRF file is the AIRS V1.0 tabulation
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
v1 = min(min(srfreq));
v2 = max(max(srfreq));

% expand interval to chunk boundaries
% v1 = 25 * floor((v1 - 5) / 25) + 5;
% v2 = 25 * ceil((v2 - 5) / 25) + 5;

