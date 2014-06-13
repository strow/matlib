function [salti, landfrac] = usgs_model_res(lat, lon, model, wgf)
% function [salti, landfrac] = usgs_model_res(lat, lon, model, wgf)
%
% model 

% Load USGS data set
%load /asl/data/usgs/world_grid_deg10_v2.mat
load(wgf)

% Make a coarse resolution version of it, based on the requested resolution:
%
% res = [res_lat res_lon ccenter]
%
% res_lat - the latitude size of the square cell
% res_lon - the longitude size of the square cell
% 
% ccenter - logical. If 0, it will assume that cell border matches the 
%           boundary -180/-90 (cell-centered grid).
%                    If 1, it will assume that cell center matches the 
%           point at -180/-90 (node-centered grid).
%

if(strcmpi(model,'ERA'))

  res = [1.5 1.5 1];


  %ilats = [0:120]'
  %ilons = [0:240];
  %mlats = +90 - ilats*1.5;
  %mlons = -180 + ilons*1.5;

elseif(strcmpi(model,'MERRA'))
 
  res = [ 1.25 1.25 1 ];

elseif(strcmpi(model,'ECMWF')) 
  res = [1./6 1./6 1];
end


% compute MODEL landfrac grid

[lf] = Image_make_map(latitude*ones(size(longitude)), ones(size(latitude))*longitude, landfrac, res);
[sa] = Image_make_map(latitude*ones(size(longitude)), ones(size(latitude))*longitude, salti, res);

ilat = floor((90-lat)/res(1)+1);
ilon = ceil((lon+180)/res(2)+1);

idx = sub2ind(size(lf),ilat, ilon);

landfrac = lf(idx);
salti = sa(idx);


end









