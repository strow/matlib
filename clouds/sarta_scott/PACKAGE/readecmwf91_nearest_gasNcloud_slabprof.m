function [head, prof, profX] = readecmwf91_nearest_gasNcloud(filename, ...
                                                             lat, lon);

% dumps in both gas profile info into prof and puts cloud info into 2 slabs,
% so that "ecmwf2cloud" and "old klayers" can be run!

% copied from /asl/matlab/gribtools/readecmwf91_nearest.m
% but now has the added functionality of including the cloud info, just like
% readecmwf91_grid_gasNcloud

%function [head, prof, profX] = readecmwf91_nearest_gasNcloud(filename, ...
%                                                                 lat, lon);
% Routine to read in a 60 or 91 level ECMWF file and return a
% RTP-like structure of profiles that are the closest grid points
% to the specified (lat,lon) locations.
%
% Input:
%    filename : (string) complete ECMWG GRIB file name
%    lat : (1 x nprof) latitudes (degrees -90 to +90)
%    lon : (1 x nprof) longitude (degrees, either 0 to 360 or -180 to 180)
%
% Output:
%    head : (RTP "head" structure of header  info) 
%    prof : (RTP "prof" structure of gas + cloud info) 
%   profX : (RTP "prof" structure of cloud   info) 
%
% Note: uses external routines: p60_ecmwf.m, p91_ecmwf.m, readgrib_inv.m,
%    readgrib_rec.m, as well as the "wgrib" program.
%

%
% Created: 19 Jan 2007  Sergio Machado
% Gas only 17 Mar 2006, Scott Hannon - re-write of old 60 level version
% Update:  23 October 2006, S.Hannon - changed from 60 level only to
%    60 or 91 level depending on file size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

xstartup
addpath /home/sergio/MATLABCODE
addpath /asl/matlab/gribtools  % for readgrib_inv.m, readgrib_rec.m
addpath /asl/matlab/aslutil    % for mktemp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the GRIB file exists
d = dir(filename);
if (length(d) ~= 1)
   disp(['did not find GRIB file=' filename]);
   return
end

tic;

% Determine 60 or 91 level based on file size.
% Note: typical file sizes are 193E+6 for 60 lev and 1160E+6 for 91 lev.
if (d.bytes < 200E+6)
   % Number of ECMWF hybrid levels
   nlev = 60;
   % Number of ECMWF latitude points
   nlat = 361;  % -90:0.50:90
   nlon = 720;  %   0:0.50:359.50
else
   % Number of ECMWF hybrid levels
   nlev = 91;
   % Number of ECMWF latitude points
   nlat = 721;  % -90:0.25:90
   nlon = 1440; %   0:0.25:359.75
end

fprintf(1,'nlev = %3i \n',nlev);

%%%%%%%%%%%%%%%%%%%
% Check lat and lon
%%%%%%%%%%%%%%%%%%%

nprof = length(lat);
if (length(lon) ~= nprof)
   disp('Error: lon and lat are different sizes!');
   return
end

% Latitude must be between -90 (south pole) to +90 (north pole)
ii = find(lat < -90 | lat > 90);
if (length(ii) > 0)
   disp('Error: latitude out of range!')
   ii
   lat(ii)
   return
end

% Note: longitude can be either 0 to 360 or -180 to 180 
ii = find(lon < -180 | lon > 360);
if (length(ii) > 0)
   disp('Error: longitude out of range!')
   ii
   lon(ii)
   return
end

% Convert any negative longitudes to positive equivalent
xlon =lon;
ii = find( lon < 0 );
xlon(ii) = 360 + lon(ii);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert (lat,lon) to fractional indices into a 2-D ECMWF grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

halfnlatp1 = round( 0.5*(nlat-1) + 1 );
nlonp1 = nlon + 1;
latres = round( (nlat-1)/180 );
lonres = round( nlon/360 );

% Convert latitude
glat = halfnlatp1 - latres*lat;
ii = find(glat < 1);
glat(ii) = 1;
ii = find(glat > nlat); % impossible except for tiny precision errors
glat(ii) = nlat;

% Convert longitude
glon = 1 + lonres*xlon;
ii = find(glon < 1);
glon(ii) = 1;
ii = find(glon >= nlonp1);  % Note: nlonp1 is 360=0 degrees
glon(ii) = 1;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the single nearest 2-D grid point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lon grid
iglon = floor( glon );
dg = glon - iglon;
ii = find( dg > 0.5);
iglon(ii) = iglon(ii) + 1;
ii = find(iglon == nlonp1);  % non-existant grid nlonp1 = grid 1 (0 deg)
iglon(ii) = 1;

% Lat grid
% Note: max(glat) corresponds to min(lat) and vice versa
iglat = floor( glat );
dg = glat - iglat;
ii = find( dg > 0.5);
iglat(ii) = iglat(ii) + 1;
clear ii dg


%%%%%%%%%%%%%%%%%%%
% 1-D ECMWF indices
%%%%%%%%%%%%%%%%%%%
% Note: in MATLAB, a 2-D (nrow, ncol) matrix is equivalent to a 1-D vector
% (nrow*ncol), with index translation 1-D index=irow + nrow*(icol-1).
%%%
% Old code
% i1D = iglat + nlat*(iglon - 1);
%%%
% C switches rows/columns compared to FORTRAN?
i1D = iglon + nlon*(iglat - 1);
clear iglat iglon

%%%
%% All needed grid points (without repeats)
%i1Dneed = unique(i1D);
%nneed = length(i1Dneed);
%
%% Create a lookup table to translate output index into index in "i1Dneed"
%i1Dtable = zeros(361*720,1);
%i1Dtable(i1Dneed) = 1:nneed;
%
%% Indices for each output profile in "i1Dneed"
%indo = i1Dtable(i1D);
%clear i1D i1Dtable
%%%

%%%
%spres=d(i1Dneed)/100;  % Divide by 100 to convert Pa to mb
%prof.spres = spres(indo);
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine record numbers for profile parameters in GRIB file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get an inventory of the GRIB file
[rec,param,level] = readgrib_inv(filename);


% Parameter "SP" surface pressure (Pa)
% Note: may be either "sfc" or "hybrid lev 1".
iparam = strcmp('SP',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find SP in GRIB inventory');
   return
end
irec_SP = ii;

% Parameter "SKT" skin temperature (K)
% Note: "sfc"
iparam = strcmp('SKT',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find SKT in GRIB inventory');
   return
end
irec_SKT = ii;

% Parameter "10U" 10 meter u wind component (m/s)
% Note: "sfc"
iparam = strcmp('10U',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find 10U in GRIB inventory');
   return
end
irec_10U = ii;

% Parameter "10V" 10 meter v wind component (m/s)
% Note: "sfc"
iparam = strcmp('10V',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find 10V in GRIB inventory');
   return
end
irec_10V = ii;

% Parameter "TCC" total cloud cover (0-1)
% Note: "sfc"
iparam = strcmp('TCC',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find TCC in GRIB inventory');
   return
end
irec_TCC = ii;

% Parameter "CI" sea ice cover (0-1)
% Note: "sfc"
iparam = strcmp('CI',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find CI in GRIB inventory');
   return
end
irec_CI = ii;

% Parameter "LSM" land/sea mask (0,1)
% Note: "sfc"; only exists in analysis files
%iparam = strcmp('LSM',param);
%ii = find(iparam == 1);
%if (length(ii) ~= 1)
%   disp('did not find LSM in GRIB inventory');
%   irec_LSM = -9999;
%else
%   irec_LSM = ii;
%end


% Parameter "T" temperature (K)
iparam_T = strcmp('T',param);
irec_T = zeros(nlev,1);

% Parameter "Q" specific humidity (kg/kg)
iparam_Q = strcmp('Q',param);
irec_Q = zeros(nlev,1);

% Parameter "O3" ozone mass mixing ratio (kg/kg)
iparam_O3 = strcmp('O3',param);
irec_O3 = zeros(nlev,1);

% Parameter "CC" cloud cover (0-1) 
iparam_CC = strcmp('CC',param);  
irec_CC = zeros(nlev,1);  
  
% Parameter "CIWC" Cloud ice water content kg/kg 
iparam_CIWC = strcmp('CIWC',param);  
irec_CIWC = zeros(nlev,1);  
  
% Parameter "CLWC" "Cloud liquid water content kg/kg 
iparam_CLWC = strcmp('CLWC',param);  
irec_CLWC = zeros(nlev,1);  

for il = 1:nlev
   strlev = ['hybrid lev ' int2str(il)];
   ilevel = strcmp(strlev,level);

   ii = find( iparam_T == 1 & ilevel == 1);
   if (length(ii) ~= 1)
      disp(['did not find T ' strlev ' in GRIB inventory']);
      return
   end
   irec_T(il) = ii;

   ii = find( iparam_Q == 1 & ilevel == 1);
   if (length(ii) ~= 1)
      disp(['did not find Q ' strlev ' in GRIB inventory']);
      return
   end
   irec_Q(il) = ii;

   ii = find( iparam_O3 == 1 & ilevel == 1);
   if (length(ii) ~= 1)
      disp(['did not find O3 ' strlev ' in GRIB inventory']);
      return
   end
   irec_O3(il) = ii;

   ii = find( iparam_CC == 1 & ilevel == 1);  
   if (length(ii) ~= 1)  
      disp(['did not find CC ' strlev ' in GRIB inventory']);  
      return  
   end  
   irec_CC(il) = ii;  
  
   ii = find( iparam_CIWC == 1 & ilevel == 1);  
   if (length(ii) ~= 1)  
      disp(['did not find CIWC ' strlev ' in GRIB inventory']);  
      return  
   end  
   irec_CIWC(il) = ii;  
  
   ii = find( iparam_CLWC == 1 & ilevel == 1);  
   if (length(ii) ~= 1)  
      disp(['did not find CLWC ' strlev ' in GRIB inventory']);  
      return  
   end  
   irec_CLWC(il) = ii;  
 
end
clear iparam iparam_T iparam_Q iparam_O3 ilevel strlev


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read surface data from the GRIB file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prof.plat = lat;
prof.plon = lon;

junk = readgrib_rec(filename,irec_SKT);
prof.stemp = junk(i1D);

junk = readgrib_rec(filename,irec_SP);
prof.spres = junk(i1D)/100; % convert Pa to hPa=mb

%if (irec_LSM > 0)
%   junk = readgrib_rec(filename,irec_LSM);
%   prof.landfrac = junk(i1D); % Note: LSM is actually a flag not a fraction
%end

% Calculate the pressure levels (using p60_ecmwf.m & p91_ecmwf.m)
prof.nlevs = nlev*ones(1,nprof);
pstr = ['prof.plevs=p' int2str(nlev) '_ecmwf( prof.spres );'];
eval(pstr);

% Assign the output header structure
head.ptype = 0;
head.pfields = 1;
head.pmin = min( prof.plevs(1,:) );
head.pmax = max( prof.plevs(nlev,:) );
head.ngas = 2;
head.glist = [1; 3];
head.gunit = [21; 21];
head.nchan = 0;
head.mwnchan = 0;

profX.plat     = prof.plat; 
profX.plon     = prof.plon; 
profX.stemp    = prof.stemp; 
profX.spres    = prof.spres; 
profX.plevs    = prof.plevs; 
%profX.landfrac = prof.landfrac;

%%%%%%%%%%%%%%%%%%%%%%%%
% Read main profile data
%%%%%%%%%%%%%%%%%%%%%%%%
prof.ptemp = zeros(nlev,nprof);
prof.gas_1 = zeros(nlev,nprof);
prof.gas_3 = zeros(nlev,nprof);

profX.ptemp = zeros(nlev,nprof);  
profX.cc    = zeros(nlev,nprof); 
profX.ciwc  = zeros(nlev,nprof); 
profX.clwc  = zeros(nlev,nprof); 

for il = 1:nlev
   junk = readgrib_rec(filename,irec_CC(il)); 
   profX.cc(il,:) = junk(i1D); 
 
   junk = readgrib_rec(filename,irec_CIWC(il)); 
   profX.ciwc(il,:) = junk(i1D); 
 
   junk = readgrib_rec(filename,irec_CLWC(il)); 
   profX.clwc(il,:) = junk(i1D); 

   junk = readgrib_rec(filename,irec_T(il));
   prof.ptemp(il,:) = junk(i1D);
   profX.ptemp(il,:) = junk(i1D); 

   junk = readgrib_rec(filename,irec_Q(il));
   prof.gas_1(il,:) = junk(i1D);
% WARNING! ECMWF water is probably specific humidity rather than mixing ratio,
% in which case this code should do: gas_1 = gas_1 / (1 - gas_1).

   junk = readgrib_rec(filename,irec_O3(il));
   prof.gas_3(il,:) = junk(i1D);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read wind data & convert
%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = readgrib_rec(filename,irec_10U);
windu = junk(i1D);
junk = readgrib_rec(filename,irec_10V);
windv = junk(i1D);

% Convert "u" and "v" wind speeds to compass direction (0-360 degrees)
% and magnitude
prof.wspeed = sqrt(windu.^2 + windv.^2);
prof.wsource = -9999*ones(1,nprof);
% WARNING!: I am not sure the direction conversion is correct
iun = find(windu < 0);
ivn = find(windv < 0);
iup = find(windu > 0);
ivp = find(windv > 0);
iu0 = find(windv == 0);
iv0 = find(windv == 0);
windu(iu0) = 1E-15;
windv(iv0) = 1E-15;
angle = atan( abs(windu)./abs(windv) )*180/pi;
% from north (ie pointing south)
ii = intersect(iu0,ivn);
prof.wsource(ii) = 0;
% from east
ii = intersect(iun,iv0);
prof.wsource(ii) = 90;
% from south
ii = intersect(iu0,ivp);
prof.wsource(ii) = 180;
% from west
ii = intersect(iup,iv0);
prof.wsource(ii) = 270;
% from 0-90 (ie pointing 180-270)
ii = intersect(iun,ivn);
prof.wsource(ii) = angle(ii);
% from 90-180 (ie pointing 270-360)
ii = intersect(iun,ivp);
prof.wsource(ii) = 180 - angle(ii);
% from 180-270 (ie pointing 0-90)
ii = intersect(iup,ivp);
prof.wsource(ii)=180 + angle(ii);
% from 270-360 (ie pointing 90-180)
ii = intersect(iup,ivn);
prof.wsource(ii) = 360 - angle(ii);
%
clear angle iun ivn iup ivp ii windu windv nneed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read remaining profile data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = readgrib_rec(filename,irec_TCC);
prof.cfrac = junk(i1D);
profX.cfrac = junk(i1D); 

junk = readgrib_rec(filename,irec_CI);
prof.udef1 = junk(i1D);

%%% end of function %%%

tnow = toc; 
fprintf(1,' took %8.6f minutes to read in ECMWF file \n',tnow/60); 
 
disp('read the ECMWF file ... parsing cloud info ... '); 

tic;
ecmwfcld2sartacld;        %%% set  the cloud profile info here 

%%%%%%%%%%%%% this is to check versus klayers_trace %%%%%%%%%%%%%%%%%%%%%

iCheck = -1;
if iCheck > 0
  check_klayers_works
  end

tnow = toc; 
fprintf(1,'TOTAL : %8.6f more minutes to process slab info \n',tnow/60); 
