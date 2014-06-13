


% Create LAT LON grid 

lat = [-90:.1:89]';
lon = [-180:.1:179];

Lat = lat*ones(size(lon));
Lon = ones(size(lat))*lon;

RTim = mattime2tai(datenum(2012,09,20))*ones(size(Lat));

prof.rtime = RTim(:);
prof.rlat = Lat(:);
prof.rlon = Lon(:);

head.ptype = 0;

pattr = set_attr('profiles','rtime','Seconds since 1993');
hattr = {};

root  = '/asl/';
[head hattr prof pattr] = rtpadd_era_data(head,hattr,prof,pattr,{'SP','SKT'})
[head hattr prof pattr] = rtpadd_ecmwf_data(head,hattr,prof,pattr,{'SP','SKT'})

% Reshape stemp
Stemp = reshape(prof.stemp,size(Lat));

figure;
scatter(flat(Lon(end-300:end,1:300)),flat(Lat(end-300:end,1:300)), 60, flat(Stemp(end-300:end,1:300)),'s','f')

% ERA is 1.5 by 1.5 - Block Cell Centered at positions: 
ilat = [0:120]';
lat_era_center = -90 + ilat*1.5;
% lat_era_bounds = -90 + [ilat-1/2 ilat+1/2]*1.5
ilon = [0:240];
lon_era_center = -180 + ilon*1.5;
% lon_era_bounds = -180 + [ilon-1/2 ilon+1/2]*1.5


% Merra:
[dat levels lats lons] = getdata_merra(datenum(2012,09,01),'ts',[],'/asl/');

[head hattr prof pattr] = rtpadd_merra(head,hattr,prof,pattr,{'SP','SKT'})




% 
dem_file='/asl/data/usgs/world_grid_deg10_v2.mat';
load(dem_file)

% compute MERRA like landfrac grid
[la15c] = Image_make_map(latitude*ones(size(longitude)), ones(size(latitude))*longitude, landfrac, [1.5 1.5 1]);

Lat_era = (-lat_era_center)*ones(size(lon_era_center));
Lon_era = ones(size(lat_era_center))*lon_era_center;

figure
scatter(flat(Lon_era(1:30,1:30)), flat(Lat_era(1:30,1:30)), 60, flat(la15c(1:30,1:30)),'s','f');



