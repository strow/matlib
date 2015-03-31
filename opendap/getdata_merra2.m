function [data levels lats lons merra_str] = getdata_merra2(time,field,level,root)
%function [data levels lats lons merra_str] = getdata_merra(time, field, level)
%
%
%  Inputs:
%       time  - matlab time
%       field - fields to be requested
%       level - data level, level field has dimensions: 1:17
%               Note, if level field is omitted or empty then 
%               all levels are retrieved
%       
%  Output
%       data  - merra data from the opendap query
%
%  Example:
%       dat = getdata_merra2(datenum(2011,3,11),'t',1);
%       simplemap(dat')
%
% 3D Fields (http://goldsmr3.sci.gsfc.nasa.gov/dods/MAI3CPASM)
%	h	geopotential height 
%	slp	sea-level pressure 
%	ps	surface pressure 
%	phis	surface geopotential 
%	o3	ozone mixing ratio 
%	qv	specific humidity 
%	ql	cloud liquid water mixing ratio 
%	qi	cloud ice mixing ratio 
%	rh	relative humidity 
%	t	air temperature 
%	u	eastward wind component 
%	v	northward wind component 
%	epv	ertel potential vorticity 
%	omega	vertical pressure velocity 
%
% 2D Fields (http://goldsmr2.sci.gsfc.nasa.gov/dods/MAT1NXSLV)
%	slp	sea level pressure 
%	ps	time averaged surface pressure 
%	u850	eastward wind at 850 hpa 
%	u500	eastward wind at 500 hpa 
%	u250	eastward wind at 250 hpa 
%	v850	northward wind at 850 hpa 
%	v500	northward wind at 500 hpa 
%	v250	northward wind at 250 hpa 
%	t850	temperature at 850 hpa 
%	t500	temperature at 500 hpa 
%	t250	temperature at 250 hpa 
%	q850	specific humidity at 850 hpa 
%	q500	specific humidity at 500 hpa 
%	q250	specific humidity at 250 hpa 
%	h1000	height at 1000 hpa 
%	h850	height at 850 hpa 
%	h500	height at 500 hpa 
%	h250	height at 250 hpa 
%	omega500	vertical pressure velocity at 500 hpa 
%	u10m	eastward wind at 10 m above displacement height 
%	u2m	eastward wind at 2 m above the displacement height 
%	u50m	eastward wind at 50 m above surface 
%	v10m	northward wind at 50 m above the displacement height 
%	v2m	northward wind at 2 m above the displacement height 
%	v50m	northward wind at 50 m above 
%	t10m	temperature at 10 m above the displacement height 
%	t2m	temperature at 2 m above the displacement height 
%	qv10m	specific humidity at 10 m above the displacement height 
%	qv2m	specific humidity at 2 m above the displacement height 
%	ts	surface skin temperature 
%	disph	displacement height 
%	troppv	pv based tropopause pressure 
%	troppt	t based tropopause pressure 
%	troppb	blended tropopause pressure 
%	tropt	tropopause temperature 
%	tropq	tropopause specific humidity 
%	cldprs	cloud-top pressure 
%	cldtmp	cloud-top temperature 
% Simple MERRA reader using internal MATLAB OpenDAP Routines

% Written by Paul Schou (Feb 2014)

if nargin < 3
  level = [];
end

switch field
  case {'slp', 'ps', 'u850', 'u500', 'u250', 'v850', 'v500', 'v250', 't850', 't500', 't250', 'q850', 'q500', 'q250', 'h1000', 'h850', 'h500', 'h250', 'omega500', 'u10m', 'u2m', 'u50m', 'v10m', 'v2m', 'v50m', 't10m', 't2m', 'qv10m', 'qv2m', 'ts', 'disph', 'troppv', 'troppt', 'troppb', 'tropt', 'tropq', 'cldprs', 'cldtmp'}
    % For the surface skin temperatures
    merra_url = 'http://goldsmr2.sci.gsfc.nasa.gov/dods/MAT1NXSLV';
    dattime = round(  (time(1) - datenum(1979,1,1,0,30,0))*24 );
    levels = 1;
  otherwise
    % For other field variables
    merra_url = 'http://goldsmr3.sci.gsfc.nasa.gov/dods/MAI3CPASM';
    dattime = round((time(1) - datenum(1979,1,1,0,0,0)) * 8);
    levels = 2;
end
    
% To open a remote dataset, use its URL:
ncid = netcdf.open ( merra_url );

% If you don't know what it contains, start by using the 'netcdf.inq' operation:
[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);

% Loop over the variables to find varid
try
  varid = netcdf.inqVarID(ncid, field);
catch
  error(['Field ''' field ''' not found in MERRA OpenDAP dataset'])
end

% Lets examine the variable:
[name,xtype,dimids,natts] = netcdf.inqVar(ncid, varid);

% Pick out the dimensions
varDims = [];
for dimid = dimids
  [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
  varDims = [varDims dimlen];
end


% Now lets get values! for var
%data = netcdf.getVar(ncid,i);
if ~isempty(level) & levels == 2
  data = netcdf.getVar(ncid, varid, [zeros(1,length(varDims)-2) level dattime], [varDims(1:end-2) 1 1]);
else
  data = netcdf.getVar(ncid, varid, [zeros(1,length(varDims)-1) dattime], [varDims(1:end-1) 1]);
end

% ATTENTION: Fill value is not consistent. Some times it is 1e15 (as
% advertised, but some times it is -9.99e8. Go figure!
data(data == -999000000) = nan;
data(data > 10^14) = nan;

if nargout > 1
  if levels > 1 % Inquire about levels if not a 2D field
    levels = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'lev'));
  end
  lats = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'lat'));
  lons = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'lon'));
  merra_str = merra_url;
end

netcdf.close(ncid);
