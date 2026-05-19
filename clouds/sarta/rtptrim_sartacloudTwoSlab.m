function [h,ha,p,pa] = rtptrim_sartacloudTwoSlab(h,ha,p,pa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iF = 0;
list = {'cprtop','cprbot','cfrac','cngwat','cpsize','ctype',...
        'cprtop2','cprbot2','cfrac2','cngwat2','cpsize2','ctype2',...
	'cfrac12',...
	'orig_ctop','orig_ctop2','sarta_lvlODice','sarta_lvlODwater'};


for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

