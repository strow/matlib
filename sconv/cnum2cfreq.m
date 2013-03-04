
function cfreq = cnum2cfreq(clist, sfile);

% NAME
%
%   cnum2cfreq -- get center frequencies from a channel list
%
% SYNOPSIS
%
%   function cfreq = cnum2cfreq(clist, sfile);
%
% INPUTS
%
%    clist  - channel number list
%    sfile  - HDF SRF tabulation data
% 
% OUTPUTS
%
%    cfreq  - channel center frequencies
% 
% DESCRIPTION
% 
%    cnum2freq takes a list of channel indices with respect to
%    a particular SRF tabulation, and returns the channel center 
%    frequencies
%
% H. Motteler, 15 Jan 02
% L. Strow, 30 Jan 04, updated srf table name
%

% the default channel list all 2378 AIRS channels
if nargin < 1
  clist = 1 : 2378;
end

% guarantee column vectors
clist = clist(:);

% the default SRF file is the V1.0 tabulation
if nargin < 2
%  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
  sfile = '/asl/data/airs/srf/srftables_031115v3.hdf';
end

% read the srf data
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

% return the associated frequencies
cfreq = freq(cind);

