function [dat levels lats lons merra_str] = getdata_merra(time, field, level)
%function [dat levels lats lons merra_str] = getdata_merra(time, field, level)
%
%
%  Inputs:
%	time  - matlab time
%	field - fields to be requested
%	level - data level, level field has dimensions: 1:17
% 	        Note, if level field is omitted or empty then 
%               all levels are retrieved
% 
%       merra_str - a string with data file information
%
%   3D Fields - 3-hr MAI3CPASM Variables: (total of 14)
%   
%   h       geopotential height
%   slp     sea-level pressure
%   ps      surface pressure
%   phis    surface geopotential
%   o3      ozone mixing ratio
%   qv      specific humidity
%   ql      cloud liquid water mixing ratio
%   qi      cloud ice mixing ratio
%   rh      relative humidity
%   t       air temperature
%   u       eastward wind component
%   v       northward wind component
%   epv     ertel potential vorticity
%   omega   vertical pressure velocity
%   
%   
%   2D Fields - Hourly MAT1NXSLV Variables: (total of 38)
%   
%   slp     sea level pressure
%   ps      time averaged surface pressure
%   u850    eastward wind at 850 hpa
%   u500    eastward wind at 500 hpa
%   u250    eastward wind at 250 hpa
%   v850    northward wind at 850 hpa
%   v500    northward wind at 500 hpa
%   v250    northward wind at 250 hpa
%   t850    temperature at 850 hpa
%   t500    temperature at 500 hpa
%   t250    temperature at 250 hpa
%   q850    specific humidity at 850 hpa
%   q500    specific humidity at 500 hpa
%   q250    specific humidity at 250 hpa
%   h1000   height at 1000 hpa
%   h850    height at 850 hpa
%   h500    height at 500 hpa
%   h250    height at 250 hpa
%   omega500 vertical pressure velocity at 500 hpa
%   u10m    eastward wind at 10 m above displacement height
%   u2m     eastward wind at 2 m above the displacement height
%   u50m    eastward wind at 50 m above surface
%   v10m    northward wind at 50 m above the displacement height
%   v2m     northward wind at 2 m above the displacement height
%   v50m    northward wind at 50 m above
%   t10m    temperature at 10 m above the displacement height
%   t2m     temperature at 2 m above the displacement height
%   qv10m   specific humidity at 10 m above the displacement height
%   qv2m    specific humidity at 2 m above the displacement height
%   ts      surface skin temperature
%   disph   displacement height
%   troppv  pv based tropopause pressure
%   troppt  t based tropopause pressure
%   troppb  blended tropopause pressure
%   tropt   tropopause temperature
%   tropq   tropopause specific humidity
%   cldprs  cloud-top pressure
%   cldtmp  cloud-top temperature
%
% See:
% http://goldsmr3.sci.gsfc.nasa.gov/dods/MAI3CPASM
% http://goldsmr2.sci.gsfc.nasa.gov/dods/MAT1NXSLV
% http://wiki.seas.harvard.edu/geos-chem/index.php/List_of_MERRA_met_fields


% Original location: /home/imbiriba/git/prod_mat/opendap/getdata_merra.m

% 3Hr files:
%	http://goldsmr3.sci.gsfc.nasa.gov/dods/MAI3CPASM/
% 1Hr files for surface skin temperature
%	http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA/MAT1NXSLV.5.2.0/YYYY/MM/MERRA300.prod.assim.tavg1_2d_slv_Nx.YYYYMMDD.hdf

  F2D_1HR = {'slp','ps','u850','u500','u250','v850','v500','v250','t850','t500','t250','q850','q500','q250','h1000','h850','h500','h250','omega500','u10m','u2m','u50m','v10m','v2m','v50m','t10m','t2m','qv10m','qv2m','ts','disph','troppv','troppt','troppb','tropt','tropq','cldprs','cldtmp'};

  F3D_3HR={'h','slp','ps','phis','o3','qv','ql','qi','rh','t','u','v','epv','omega'};

  CLD_3D_3HR={'rh','qlls','qils','qlan','qian','qccu','cfls','cfan','cfcu'};


  switch field
    %case {'ts','u2m','v2m'}
    case F2D_1HR
      % Call for the surface skin temperature 
      dattime_1h = round(  (time(1) - datenum(1979,1,1,0,30,0))*24 );
      filename = ['/asl/data/merra/' datestr(time(1),'yyyy/mm/dd') '/MAT1NXSLV_' field '_' datestr(round(time(1)*24)/24,'yyyymmdd-HHMMSS') '.mat'];

    case F3D_3HR
      % Call for other field variables
      dattime = round((time(1) - datenum(1979,1,1,0,0,0)) * 8);
      filename = ['/asl/data/merra/' datestr(time(1),'yyyy/mm/dd') '/MAI3CPASM_' field '_' datestr(round(time(1)*8)/8,'yyyymmdd-HHMMSS') '.mat'];

    case CLD_3D_3HR
      dattime = round((time(1) - datenum(1979,1,1,0,0,0)) * 8);
      filename = ['/asl/data/merra/' datestr(time(1),'yyyy/mm/dd') '/MAT3CPCLD_' field '_' datestr(round(time(1)*8)/8,'yyyymmdd-HHMMSS') '.mat'];
  
    otherwise
      error('getdata_merra: Unknown requested field');
  end

  try
    mkdirs(dirname(filename));
  catch
    disp(['Cannot create ' dirname(filename)]);
    return
  end


  merra_str = filename;
  if exist(filename,'file')
      levels = [];
      load(filename)
      if nargin == 2 || length(level) == 0
	return
      end
      dat = dat(:,:,level);
      return
  end


  try
  switch field
    %case {'ts','u2m','v2m'}
    case F2D_1HR
      url = 'http://goldsmr2.sci.gsfc.nasa.gov/dods/MAT1NXSLV';
      fld = [field '[' num2str(dattime_1h) ':1:' num2str(dattime_1h) '][0:360][0:539]'];

      [dat x lats lons] = getdata_opendap(url,fld);

      % NOTE: The returned time can be converted to matlab time by doing this:
      % mtimes = times+365;  % i.e. it is one year 
      % mtimes = x + 365;
      levels=[];
      save(filename,'dat','lats','lons');
  
    % "3D" that are actually 2D
    case {'slp','ps','phis'}
      url = 'http://goldsmr3.sci.gsfc.nasa.gov/dods/MAI3CPASM';
      fld = [field '[' num2str(dattime) ':1:' num2str(dattime) '][0:1:143][0:1:287]'];
      [dat x lats lons] = getdata_opendap(url, fld);

      levels = [];
      save(filename,'dat','lats','lons');

    % Remaining 3D fields
    case F3D_3HR
      url = 'http://goldsmr3.sci.gsfc.nasa.gov/dods/MAI3CPASM';
      fld = [field '[' num2str(dattime) ':1:' num2str(dattime) '][0:1:41][0:1:143][0:1:287]'];
      [dat x levels lats lons] = getdata_opendap(url, fld);

      save(filename,'dat','levels','lats','lons');

      if nargin == 2 || length(level) == 0
	return
      end
      dat = dat(:,:,level);

    case CLD_3D_3HR
      url = 'http://goldsmr3.sci.gsfc.nasa.gov/dods/MAT3CPCLD';
      fld = [field '[' num2str(dattime) ':1:' num2str(dattime) '][0:1:41][0:1:143][0:1:287]'];

      [dat x levels lats lons] = getdata_opendap(url, fld);

      save(filename,'dat','levels','lats','lons');

      if nargin == 2 || length(level) == 0
	return
      end
      dat = dat(:,:,level);

    otherwise 
      error(['You cant really be here!! Bad error']);
  end

  catch thiserr
    estr={};
    estr{1}=('ERROR in GETDATA_MERRA while fetching OpenDap data:');
    estr{2}='.';
    estr{3}=(thiserr.message);
    for ii=1:numel(thiserr.stack)
      estr{end+1}=(['.      file: ' thiserr.stack(ii).file]);
      estr{end+1}=(['.      name: ' thiserr.stack(ii).name]);
      estr{end+1}=(['.      line: ' num2str(thiserr.stack(ii).line)]);
      estr{end+1}='.';
    end  
    estr{end+1}=('Requested URL:');
    estr{end+1}=([url '.dods?' fld]);
    error(char(estr)'); 
  end 
  

end
