
function [flist, dout] = modskip(clist, din)

% NAME
%
%   modskip -- add NaNs between modules for "pen lift"
%
% SYNOPSIS
%
%   [flist, dout] = modskip(clist, din)
%   
% INPUTS
%
%    clist  - m-vector of channel numbers
%    din    - m x n array of data
%
% OUTPUTS
%
%    flist  - channel freq's, with NaN's at module boundaries
%    dout   - din, with NaN's added to rows at module boundaries
%
% DESCRIPTION
% 
%    modskip takes an m-vector of channel numbers and an m x n data
%    array and returns an (m+k) x n data array, with k extra rows
%    of NaNs are added between module boundaries.  Module boundaries
%    are read from a matlab file "clist.mat".
%

clist = clist(:);
nchan = length(clist);
[m,ncol] = size(din);
if m ~= nchan
  error('number of rows of din must match length of clist');
end

% make sure channel numbers are in increasing order
[clist, p1] = sort(clist);

% keep data in channel order
din = din(p1, :);

% get up to date channel and module info:
%
%    cfreq    - channel center frequency
%    cnum     - "official" channel number
%    cmod     - channel module string names (one per channel)
%    modlist  - list of all module string names
%
% the channel set is sorted by module frequency bands, and then 
% by frequency within modules
%
% (this needs to be in a standard place)

load /asl/data/airs/srf/clist

% get indices of the channels specified in clist
cind = interp1(cnum, 1:length(cnum), clist, 'nearest');

% get frequencys and modules for the requested channels
cfreq2 = cfreq(cind);
cmod2 = cmod(cind);

% initialize the index into the output arrays
dind = zeros(nchan, 1);
dind(1) = 1;
k = 1;

% build the index, add a skip at each module boundary
for i = 2 : nchan
  k = k + 1;
  if ~isequal(cmod2(i-1), cmod2(i))
    k = k + 1;
  end
  dind(i) = k;
end

% initialize the output arrays and copy the input data to them
flist = zeros(k, 1);
dout = zeros(k, ncol);

flist = flist .* NaN;
flist(dind) = cfreq2;

dout = dout .* NaN;
dout(dind, :) = din;

