function [h ha p pa] = rtpadd_usgs_10dem(h,ha,p,pa)
% function [h,ha,p,pa] = rtpadd_usgs_10dem(h,ha,p,pa)
%
% Add surface altitude and land fraction based on the usgs_10dem.m data.
%
% Breno Imbiriba - 2013.03.14


  [salti landfrac] = usgs_deg10_dem(p.rlat, p.rlon);

  p.salti=salti;
  p.landfrac=landfrac;

  ha=set_attr(ha,'topo','usgs_deg10_dem');

end

