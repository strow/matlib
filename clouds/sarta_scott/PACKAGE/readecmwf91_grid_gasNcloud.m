function [head, prof, profX] = readecmwf91_grid_gasNcloud(filename, ...
                                  minlat, maxlat, minlon, maxlon, gridsize);

% function [head, prof, profX] = readecmwf91_grid_gasNcloud(filename, ...
%    minlat, maxlat, minlon, maxlon, gridsize);
%
% Routine to read in ECMWF file and return a RTP-like structure of
% profiles for the rectangular area specified by the min/max lat/lon
% locations.  gridsize may be used to thin the grid.  The min/max
% boundaries and gridsize must be choosen so that the output starts
% and ends on the specified boundaries.
%
% Input:
%    filename : (string) complete ECMWG GRIB file name
%    minlat   : (scaler) mininum latitudes (-90 to +90 degrees)
%    maxlat   : (scaler) maximum latitudes (-90 to +90 degrees)
%    minlon   : (scaler) mininum longitude (0 to 359.75 degrees)
%    maxlon   : (scaler) maximum longitude (0 to 359.75 degrees)
%    gridsize : (scaler) grid size, ie spacing of points (degrees)
%
% Output:
%    head : (RTP "head" structure of header  info)
%    prof : (RTP "prof" structure of gas + cloud info)
%   profX : (RTP "prof" structure of cloud   info)
%
% Note: uses external routines p_ecmwf.m and svp.m
%

% Created: 17 Mar 2006, Scott Hannon - re-write of old 60 level code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

addpath /asl/matlab/gribtools  % for readgrib_inv.m, readgrib_rec.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

% Number of ECMWF hybrid levels
nlev = 91;

% Number of ECMWF latitude points
nlat = 721;  % -90:0.25:90
nlon = 1440; %   0:0.25:359.75

halfnlatp1 = round( 0.5*(nlat-1) + 1 );
nlonp1 = nlon + 1;
latres = round( (nlat -1)/180 );
lonres = round( nlon/360 );
maxlongrid = 0.0001 + 360 - 1/lonres;
if (latres ~= lonres)
   disp('Error: latitude resolution differs from longitude resolution')
end


%%%%%%%%%%%%%%%%%%%
% Check lat and lon
%%%%%%%%%%%%%%%%%%%
% Latitude
if (minlat < -90 | minlat > 90)
   disp('Error: minlat outside allowed range -90 to +90')
   return
end
if (maxlat < -90 | maxlat > 90)
   disp('Error: minlat outside allowed range -90 to +90')
   return
end
if (minlat > maxlat)
   disp('Error: minlat > maxlat')
   return
end

% Longitude
if (minlon < 0 | minlon > maxlongrid)
   disp(['Error: minlon outside allowed range 0 to ' num2str(maxlongrid)])
   return
end
if (maxlon < 0 | maxlon > maxlongrid)
   disp(['Error: maxlon outside allowed range 0 to ' num2str(maxlongrid)])
   return
end
if (minlon > maxlon)
   disp('Error: minlon > maxlon')
   return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert min/max lat/lon into grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Latitude
% Note since grid run 90 to -90, imaxglat corresponds to min (not max) lat
imaxglat = round( halfnlatp1 - latres*minlat );
iminglat = round( halfnlatp1 - latres*maxlat );

% Longitude
iminglon = round( 1 + lonres*minlon );
imaxglon = round( 1 + lonres*maxlon );


%%%%%%%%%%%%%%%%
% Check gridsize
%%%%%%%%%%%%%%%%
igridstep = round( latres*gridsize );

% Latitude grid
xglatsize = 1 + (imaxglat - iminglat)/igridstep;
iglatsize = round( xglatsize );
if ( abs(iglatsize - xglatsize) > 0.001 )
   disp('Error: gridsize not compatible with min/max lat');
   return
end

% Longitde grid
xglonsize = 1 + (imaxglon - iminglon)/igridstep;
iglonsize = round( xglonsize );
if ( abs(iglonsize - xglonsize) > 0.001 )
   disp('Error: gridsize not compatible with min/max lon');
   return
end
clear xglonsize xglatsize


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create matrices with lat/lon of complete grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dummylat = zeros(nlat,nlon);
dummylon = zeros(nlat,nlon);
res = 1/latres;
junklat = 90:-res:-90;
junklon = 0:res:(360-res);
junkones=ones(1,nlon);
for ii=1:nlat
   dummylat(ii,:) = junklat(ii).*junkones;
   dummylon(ii,:) = junklon;
end
clear junklat junklon junkones ii


%%%%%%%%%%%%%%%%%%%%%%%
% Determine 1-D indices
%%%%%%%%%%%%%%%%%%%%%%%
% Note: in MATLAB, a 2-D (nrow, ncol) matrix is equivalnt to a 1-D vector
% (nrow*ncol), with index translation 1-D index=irow + nrow*(icol-1).

%%%
% Old code
% rowterm = round( iminglat:igridstep:imaxglat );
% colterm = round( 361 * ((iminglon:igridstep:imaxglon) - 1) );
% index1d = rowterm'*ones(1,iglonsize) + ones(iglatsize,1)*colterm;
%%%
% C switches row/columns compared to FORTRAN?
colterm = round( iminglon:igridstep:imaxglon );
rowterm = round( nlon * ((iminglat:igridstep:imaxglat) - 1) )';
index1d = ones(iglatsize,1)*colterm + rowterm*ones(1,iglonsize,1);
dummylat = dummylat';
dummylon = dummylon';

nprof = round( iglonsize * iglatsize);
i1D = reshape(index1d, 1, nprof);
clear colterm rowterm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine record numbers for profile parameters in GRIB file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the GRIB file exists
d = dir(filename);
if (length(d) ~= 1)
   disp(['did not find GRIB file=' filename]);
   return
end


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
iparam = strcmp('LSM',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find LSM in GRIB inventory');
   irec_LSM = -9999;
else
   irec_LSM = ii;
end


% Parameter "Z" geopotential (m^2/s^2)
% Note: only exists in analysis files where it might appear twice
% as "sfc" and/or "hybrid lev 1".
iparam = strcmp('Z',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   if (length(ii) == 2)
      irec_Z = ii(1);
   else
      disp('did not find Z in GRIB inventory');
      irec_Z = -9999;
   end
else
   irec_Z = ii;
end

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
clear        iparam_CC iparam_CIWC iparam_CLWC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read surface data from the GRIB file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prof.plat = dummylat(i1D);
prof.plon = dummylon(i1D);

junk = readgrib_rec(filename,irec_SKT);
prof.stemp = junk(i1D);

junk = readgrib_rec(filename,irec_SP);
prof.spres = junk(i1D)/100; % convert Pa to hPa=mb

if (irec_LSM > 0)
   junk = readgrib_rec(filename,irec_LSM);
   prof.landfrac = junk(i1D); % Note: LSM is actually a flag not a fraction
end
if (irec_Z > 0)
   junk = readgrib_rec(filename,irec_Z);
   prof.salti = junk(i1D)/9.8; % convert geopotential to altitude
end


% Calculate the pressure levels
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
profX.landfrac = prof.landfrac;
profX.plevs    = prof.plevs;

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
ecmwfcld2sartacld;        %%% set  the cloud profile info here 