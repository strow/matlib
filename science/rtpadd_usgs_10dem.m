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
%
% Breno Imbiriba - 2013.03.14

  if(~exist('root','var'))
    wgf = '/asl/data/usgs/world_grid_deg10.mat'
  else
    wgf = [root '/data/usgs/world_grid_deg10.mat'];
  end

  [salti landfrac] = usgs_deg10_dem(p.rlat, p.rlon,wgf);

  p.salti=single(salti);
  p.landfrac=single(landfrac);

  ha=set_attr(ha,'topo','usgs_deg10_dem');

end

