function [rout, fout] = sconv2_vd(rin, fin, clist, vdsfile, coff);

% function [rout, fout] = sconv2_vd(rin, fin, clist, vdsfile, coff);
%
% Do convolutions directly from tabulated SRFs.  This "_vd" version
% of sconv2 uses Paul van Delst's netcdf format of 0.0025 cm^-1
% gridded SRFs rather than UMBC's HDF format.
%
% Input:
%    rin    - [m x n] radiance/transmittance data to be convolved
%    fin    - [m x 1] frequencies points corresponding to rin
%    clist  - [k x 1] desired channel IDs
%    vdsfile- (string) filename of van Delst's netcdf format SRF data
%    coff   - optional [k x 1] channel center offsets
% 
% Output:
%    rout   - [k x n] convolved data
%    fout   - [k x 1] approximate channel center frequencies (including
%                effects of coff)

% Created (sconv2.m) : 15 Jan 2002, H.Motteler
% Created (sconv2_vd.m): 04 Jan 2005, S.Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlab/h4tools

% guarantee input vectors are columns
fin = fin(:);
clist = clist(:);


% Check clist is unique (no repeats)
nchan = length(clist);
junk = unique( clist );
if (length(junk) ~= nchan)
   error('each channel ID in clist must be unique; no repeats allowed')
end


% check rin dimensions
[m,n] = size(rin);
if m ~= length(fin)
  error('number of rows of rin do not match length of fin')
end


% Check SRF file
if nargin < 4
   error('insufficient input arguments were specified')
   return
end
d = dir(vdsfile);
if (length(d) ~= 1)
   error(['unable to read vdsfile=' vdsfile])
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


% read the SRF data
[dat, atr, dim] = sdload(vdsfile);
clear atr dim
chanid = double( dat.channel_list );


% Check SRF file contains all channels specified by clist
[junk, ia] = setdiff( clist, chanid );
if (length(junk) > 0)
   error(['clist contains IDs not found in vdsfile: ' int2str(junk)]);
end


% initialize output array
rout = zeros(nchan, n);
fout = zeros(nchan, 1);

% loop on SRFs
for j = 1 : nchan

  ic = clist(j);            % current channel ID
  ii = find(chanid == ic);  % index of channel in SRF freq arrays

  eval(['srfval = dat.channel_' int2str(ic) '_response'';']);
  npts = length(srfval);

  % get the frequency span of the current SRF
  v1 = round( dat.begin_frequency(ii)*400 )/400;
  v2 = round( dat.end_frequency(ii)*400 )/400;
  %
  % Apply channel offset (if any)
  if (coff(j) ~= 0)
     v1 = v1 + coff(j);
     v2 = v2 + coff(j);
  end
  %
  srfreq = v1:0.0025:v2;
  if (length(srfreq) ~= npts)
     npts
     error(['channel ' int2str(ic) '  length of srfreq & srfval differ']);
  end

  % if the SRF is entirely outside fin, skip this channel
  if v2 <= f1 | f2 <= v1
    fprintf(1,'sconv2_vd(): WARNING -- SRF %d outside of input range\n',ic);
    continue
  end

  % if the SRF overlaps fin, lop it off to fit; 
  % give a warning message if we lop more than dv
  if v1 < f1
    if v1 < f1 - dv
      fprintf(1,'sconv2vd(): WARNING -- truncating LHS of SRF %d\n',ic);
    end
    v1 = f1;
  end
  if f2 < v2
    if f2 + dv < v2
      fprintf(1,'sconv2_vd(): WARNING -- truncating RHS of SRF %d\n',ic);
    end
    v2 = f2;
  end

  % find the indices of the current SRF in fin 
  v1ind = ceil((v1-f1)/dv) + 1;
  v2ind = floor((v2-f1)/dv) + 1;

  % interpolate the SRF to a subinterval of the fin grid
  s1 = interp1(srfreq, srfval, fin(v1ind:v2ind), 'linear');

  % normalize the SRF
  s1sum = sum(s1);
  s1cumsum = cumsum(s1);
  s1 = s1 ./ s1sum;
  
  % Find area weighted channel mid-point
  fsrf = fin(v1ind:v2ind);
  frac = s1cumsum./s1sum;
  ilo = max( find( frac <= 0.5 ) );
  ihi = ilo + 1;
  fout(j) = ( (fsrf(ihi) - fsrf(ilo))./(frac(ihi)-frac(ilo)) ) .* ...
     (0.5 - frac(ilo)) + fsrf(ilo);
  
  % apply the SRF
  rout(j, :) = s1' * rin(v1ind:v2ind, :);
end

%%% end of function %%%
