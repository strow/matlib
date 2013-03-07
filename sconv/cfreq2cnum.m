
function clist = cfreq2cnum(cfreq, sfile);

% NAME
%
%   cfreq2cnum -- get a channel list from center frequencies
% 
% SYNOPSIS
%
%   function clist = cfreq2cnum(cfreq, sfile);
%
% INPUTS
%
%   cfreq  - channel frequencies
%   sfile  - HDF SRF tabulation data
% 
% OUTPUTS
%
%   clist   - channel numners
% 
% DESCRIPTION
% 
%   cfreq2num takes a list of channel frequencies and returns a
%   channel set, as indices with respect to a particular SRF tabulation.  
%   For each frequency in cfreq, the closest matching frequency in the
%   SRF tabulation is found, and the associated channel number returned.
%   
% H. Motteler, 15 Jan 02
% L. Strow, 30 Jan 04, updated srf table name
%


% guarantee column vectors
cfreq = cfreq(:);

% set a default SRF file
if nargin < 2
%  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
%  sfile = '/asl/data/airs/srf/srftables_m135_fringes_nov02.hdf';
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

[nchan, nspts] = size(srfval);

% match supplied and tabulated frequencies
cind = interp1(freq, 1:nchan, cfreq, 'nearest');

% return associated channel numbers
clist = chanid(cind);

