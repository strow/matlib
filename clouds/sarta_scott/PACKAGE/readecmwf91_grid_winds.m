function [head, profX] = readecmwf91_grid_winds(filename, minlat, maxlat, ...
   minlon, maxlon, gridsize);

% function [head, profX] = readecmwf91_grid_winds(filename, minlat, maxlat,...
%    minlon, maxlon, gridsize);
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
%    head : (RTP "head" structure of header info)
%    profX : (RTP "profX" structure of profile info)
%
% Note: uses external routines p_ecmwf.m and svp.m
%

% Created: 07 Jan 2007, Sergio Machado -- readecmwf91_grid.m for cloud fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

addpath /asl/matlab/gribtools  % for readgrib_inv.m, readgrib_rec.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Parameter "2T" 2m temperature
% Note: may be either "sfc" or "hybrid lev 1". 
iparam = strcmp('2T',param); 
ii = find(iparam == 1); 
if (length(ii) ~= 1) 
   disp('did not find 2T in GRIB inventory'); 
   return 
end 
irec_2T = ii; 

% Parameter "10U" surface U speed
% Note: may be either "sfc" or "hybrid lev 1". 
iparam = strcmp('10U',param); 
ii = find(iparam == 1); 
if (length(ii) ~= 1) 
   disp('did not find 10U in GRIB inventory'); 
   return 
end 
irec_10U = ii; 

% Parameter "10V" surface V speed
% Note: may be either "sfc" or "hybrid lev 1". 
iparam = strcmp('10V',param); 
ii = find(iparam == 1); 
if (length(ii) ~= 1) 
   disp('did not find 10V in GRIB inventory'); 
   return 
end 
irec_10V = ii; 

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

% Parameter "TCC" total cloud cover (0-1)
% Note: "sfc"
iparam = strcmp('TCC',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find TCC in GRIB inventory');
   return
end
irec_TCC = ii;

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
 
% Parameter "T" temperature (K) 
iparam_T = strcmp('T',param); 
irec_T = zeros(nlev,1); 

% Parameter "CC" cloud cover (0-1)
iparam_CC = strcmp('CC',param); 
irec_CC = zeros(nlev,1); 
 
% Parameter "U" speed
iparam_U = strcmp('U',param); 
irec_U = zeros(nlev,1); 
 
% Parameter "V" speed
iparam_V = strcmp('V',param); 
irec_V = zeros(nlev,1); 
 
for il = 1:nlev 
   strlev = ['hybrid lev ' int2str(il)]; 
   ilevel = strcmp(strlev,level); 

   ii = find( iparam_T == 1 & ilevel == 1); 
   if (length(ii) ~= 1) 
      disp(['did not find T ' strlev ' in GRIB inventory']); 
      return 
   end 
   irec_T(il) = ii; 
 
   ii = find( iparam_CC == 1 & ilevel == 1); 
   if (length(ii) ~= 1) 
      disp(['did not find CC ' strlev ' in GRIB inventory']); 
      return 
   end 
   irec_CC(il) = ii; 
 
   ii = find( iparam_U == 1 & ilevel == 1); 
   if (length(ii) ~= 1) 
      disp(['did not find U ' strlev ' in GRIB inventory']); 
      return 
   end 
   irec_U(il) = ii; 
 
   ii = find( iparam_V == 1 & ilevel == 1); 
   if (length(ii) ~= 1) 
      disp(['did not find V ' strlev ' in GRIB inventory']); 
      return 
   end 
   irec_V(il) = ii; 
 
end 
clear iparam iparam_CC iparam_U iparam_V ilevel strlev

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read surface data from the GRIB file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

profX.plat = dummylat(i1D);
profX.plon = dummylon(i1D);

junk = readgrib_rec(filename,irec_2T);
profX.T2 = junk(i1D);

junk = readgrib_rec(filename,irec_10U);
profX.U10 = junk(i1D);

junk = readgrib_rec(filename,irec_10V);
profX.V10 = junk(i1D);

junk = readgrib_rec(filename,irec_SKT);
profX.stemp = junk(i1D);

junk = readgrib_rec(filename,irec_SP);
profX.spres = junk(i1D)/100; % convert Pa to hPa=mb

if (irec_LSM > 0)
   junk = readgrib_rec(filename,irec_LSM);
    profX.landfrac = junk(i1D); % Note: LSM is actually a flag not a fraction
end

% Calculate the pressure levels
profX.nlevs = nlev*ones(1,nprof);
pstr = ['profX.plevs=p' int2str(nlev) '_ecmwf( profX.spres );'];
eval(pstr);

% Assign the output header structure
head.ptype = 0;
head.pfields = 1;
head.pmin = min( profX.plevs(1,:) );
head.pmax = max( profX.plevs(nlev,:) );
head.ncld = 3;               %%% 3 cloud parameters being output
head.gcld = [1;2;3];         %%% CC, U, V
head.nchan = 0;
head.mwnchan = 0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Read main profile data
%%%%%%%%%%%%%%%%%%%%%%%%
profX.ptemp = zeros(nlev,nprof); 
profX.cc = zeros(nlev,nprof);
profX.uspd = zeros(nlev,nprof);
profX.vspd = zeros(nlev,nprof);
for il = 1:nlev
   junk = readgrib_rec(filename,irec_T(il)); 
   profX.ptemp(il,:) = junk(i1D); 
   
   junk = readgrib_rec(filename,irec_CC(il));
   profX.cc(il,:) = junk(i1D);

   junk = readgrib_rec(filename,irec_U(il));
   profX.uspd(il,:) = junk(i1D);

   junk = readgrib_rec(filename,irec_V(il));
   profX.vspd(il,:) = junk(i1D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read remaining profile data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = readgrib_rec(filename,irec_TCC);
profX.cfrac = junk(i1D);

%%% end of function %%%
