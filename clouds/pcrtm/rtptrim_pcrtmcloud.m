function [h,ha,p,pa] = rtptrim_pcrtmcloud(h,ha,p,pa)

[h,ha,p,pa] = rtptrim_sartacloud(h,ha,p,pa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iF = 0;
clear list
list = {'ncol','overlap','pcrtm_co2_used','pcrtm_iceODX','pcrtm_waterODX',...
        'pcrtm_lvlODice','pcrtm_lvlODwater','pcrtm_iceODX','pcrtm_waterODX',...
        'pcrtm_iceOD','pcrtm_iceDME','pcrtm_iceCTOP','pcrtm_waterOD','pcrtm_waterDME','pcrtm_waterCTOP'....
	'rcalc_std','rad_allsky_std'};
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end

list = {'sarta_xclear','sarta_rclearcalc','sarta_lvlODice','sarta_lvlODwater',...
	'sarta_clr_co2_used','sarta_xclr_co2_used','sarta_cld_co2_used'};
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end
