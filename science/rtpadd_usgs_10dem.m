function [h ha p pa] = rtpadd_usgs_10dem(h,ha,p,pa,root)
% function [h,ha,p,pa] = rtpadd_usgs_10dem(h,ha,p,pa,root)
%
%   h,ha,p,pa - rtp structure
% 
%   Load the USGS matlab file: /asl/data/usgs/world_grid_deg10.mat
%
%   root - optional: root directory (usually /asl) of where to find 
%                    the USGS data file. 
%
% Add surface altitude and land fraction based on the usgs_10dem.m data.
% Invalid Lats and Lons will be ignored and marked as -9999;
%
% Breno Imbiriba - 2013.03.14

  if(~exist('root','var'))
    wgf = '/asl/data/usgs/world_grid_deg10.mat'
  else
    wgf = [root '/data/usgs/world_grid_deg10.mat'];
  end

  % If there's any bad GEO data, replace it by the (0,0) so not to crash
  % usgs_deg10_dem.m. 
  ibad_geo = find(abs(p.rlat)>90 | p.rlon<-180 | p.rlon>360);
  bad_rlat = p.rlat(ibad_geo);
  bad_rlon = p.rlon(ibad_geo);
  p.rlat(ibad_geo) = 0; 
  p.rlon(ibad_geo) = 0;

  % Call main geo routine
  [salti landfrac] = usgs_deg10_dem(p.rlat, p.rlon,wgf);

  p.salti=single(salti);
  p.landfrac=single(landfrac);

  % Replace bad points by -9999
  p.rlat(ibad_geo) = bad_rlat;
  p.rlon(ibad_geo) = bad_rlon;
  p.salti(ibad_geo)    = -9999;
  p.landfrac(ibad_geo) = -9999;


  ha=set_attr(ha,'topo','usgs_deg10_dem');

end

