function [head, hattr, prof, pattr, current_ename] = rtpadd_grib_data(current_ename, sourcename, head, hattr, prof, pattr);

% function [head, hattr, prof, pattr] = rtpadd_grib_data(sourcename, head, hattr, prof, pattr);
%
% Routine to read in a grib file and return a RTP-like structure of
% profiles that are the closest grid points to the specified (lat,lon)
% locations, or interpolated from the grib data (default)
%
% Input:
%    sourcename : (string) complete ECMWG GRIB file name
%                to automatically pick files use either 'ECMWF' or 'ERA'
%    head      : rtp header structure
%    hattr     : header attributes
%    prof.       profile structure with the following fields
%        rlat  : (1 x nprof) latitudes (degrees -90 to +90)
%        rlon  : (1 x nprof) longitude (degrees, either 0 to 360 or -180 to 180)
%        rtime : (1 x nprof) observation time in seconds
%    pattr     : profile attributes, note: rtime must be specified
%
% Output:
%    head : (RTP "head" structure of header info)
%    hattr: header attributes
%    prof : (RTP "prof" structure of profile info)
%    pattr: profile attributes
%
% Below mfile usage needs to be updated
% Note: uses external routines: p60_ecmwf.m, p91_ecmwf.m, readgrib_inv_data.m,
%    readgrib_offset_data.m, as well as the "wgrib" program.
%

% Created: 17 Mar 2006, Scott Hannon - re-write of old 60 level version
% Rewrite:  4 May 2011, Paul Schou - switched to matlab binary reader
% Update : 17 Jun 2011, Paul Schou - added grib 2 capabilities
% Update: 27 Jun 2011, S.Hannon - add isfield test for head.pfields
% Update: 05 Jan 2012, L. Strow - bitor argument translated to uint32
% June 20, 2014.  Complete re-write by L. Strow, now using netcdf version of 
%     grib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Taken from fill_ecmwf.m in ~strow/Matlab/Gribools

% Location of grib files (add root var later)
fhdr = fullfile('/asl/data/',lower(sourcename));

% Now input parameter
% current_ename = '';  % This should be placed outside a rtp file loop
% Load new rtp file here

% Hardcoded for CrIS and IASI (iasi2mattime)
mtime = iasi2mattime(prof.rtime);

% Get a cell array of ecmwf grib files for each time
enames = get_enames(mtime);

% Find the unique grib files and indices that go with them
[u_enames, ia, ic] = unique(enames);
n = length(u_enames);

% Loop over unique grib file names
for i = 1:n
% Build file name from parts
   fne = ['UAD' u_enames{i} '001'];
   e_mth_year = datestr(mtime(ia(i)),'yyyymm');
   fn = fullfile(fhdr,e_mth_year(1:4),e_mth_year(5:6),fne);
% Actually read grib1, grib2 .nc files
   fn_s = [fn '-1.nc'];
   fn_h = [fn '-2.nc'];
% Do the netcdf files exist?
   if exist(fn_s) == 0 | exist(fn_h) == 0 
      disp(['Netcdf grib files missing for root: ' fn])
      break % Go to next partition
   end
% If the filename has changed, re-load F   
   if ~strcmp(current_ename,fn) 
      F = grib_interpolate(fn_s,fn_h);
      current_ename = fn;
   end   
% Fill rtp fields
   k = find( ic == i );  % indices of first partition (of n total)
% Assume rtp lat/lon are +-180??  Need to be 0-360 for grib interpolation
   rlat = prof.rlat(k);
   rlon = prof.rlon(k);
   rlon(rlon<0) = rlon(rlon<0) + 360;

   prof.sstk(k)    = F.sst.ig(rlat,rlon);
   prof.spres(k)   = F.sp.ig(rlat,rlon);
   prof.stemp(k)   = F.skt.ig(rlat,rlon);
   wind_v          = F.v10.ig(rlat,rlon);
   wind_u          = F.u10.ig(rlat,rlon);
   prof.wspeed(k)  = sqrt(wind_u.^2 + wind_v).^2;
   prof.wsource(k) = mod(atan2(single(wind_u), single(wind_v)) * 180/pi,360);
   prof.cfrac(k)   = F.tcc.ig(rlat,rlon);
   ci_udef = 1;
   prof.udef(ci_udef,k) = F.ci.ig(rlat,rlon);
% F.tcwv.ig  % Total column water?  Use this instead of ours?
% F.msl.ig   % Not in rtp for now

% Hybrid parameters
% levid = 1 is top of atmosphere
% b are the sortedd level IDs   
   [b,j]=sort(F.levid);
   for l=1:length(F.levid)
      prof.ptemp(l,k) = F.t(j(l)).ig(rlat,rlon);
      prof.gas_1(l,k) = F.q(j(l)).ig(rlat,rlon);
      prof.gas_3(l,k) = F.o3(j(l)).ig(rlat,rlon);
      prof.cc(l,k)    = F.cc(j(l)).ig(rlat,rlon);
      prof.clwc(l,k)  = F.clwc(j(l)).ig(rlat,rlon);
      prof.ciwc(l,k)  = F.ciwc(j(l)).ig(rlat,rlon);
   end
% Below hard-coded for now, put all level calcs into one file
% Use max(levid) to tell routine which model   
% Only want pressure levels in grib file, in order
   xtemp = p137_ecmwf(prof.spres(k));  % all 137 pressure levels
   prof.plevs(:,k) = xtemp(b,:);  % subset to ones in grib file
end

% Header info
head.ptype = 0;
% Do I need stuff below commented out?
% if (isfield(head,'pfields'))
%    head.pfields = bitor(uint32(head.pfields), 1);
% else
%    head.pfields = 1;
% end
head.ngas = 2;
head.glist = [1; 3];
head.gunit = [21; 21];
head.pmin = min( prof.plevs(1,:) );
head.pmax = max( prof.plevs(end,:) );
% Setting attributes needs work...
% pattr = set_attr(pattr,'profiles','ECMWF','profiles');

% I think this is needed to avoid negatives in SARTA?
min_H2O_gg = 3.1E-7;  % 0.5 pppm
min_O3_gg = 1.6E-8;   % 0.01 ppm
% Find/replace bad mixing ratios
if isfield(prof,'gas_1')
  ibad = find(prof.gas_1 <= 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.gas_1(ibad) = min_H2O_gg;
%    say(['Replaced ' int2str(nbad) ' negative/zero H2O mixing ratios'])
  end
end
%
if isfield(prof,'gas_3')
  ibad = find(prof.gas_3 <= 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.gas_3(ibad) = min_O3_gg;
%    say(['Replaced ' int2str(nbad) ' negative/zero O3 mixing ratios'])
  end
end
%  fix any cloud frac
if isfield(prof,'cfrac')
  ibad = find(prof.cfrac > 1);
  nbad = length(ibad);
  if (nbad > 0)
    prof.cfrac(ibad) = 1;
%    say(['Replaced ' int2str(nbad) ' CFRAC > 1 fields'])
  end
end

