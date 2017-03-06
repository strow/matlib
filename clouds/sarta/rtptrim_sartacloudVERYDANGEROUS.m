function [h,ha,p,pa] = rtptrim_sartacloudVERYDANGEROUS(h,ha,p,pa)

%% this function is similar to rtptrim_sartacloud.m, and does more trimming
%% of the cloud fioleds that are produced

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iF = 0;	

clear list
list = {'sarta_lvlODice','sarta_lvlODwater','sarta_xclear','orig_ctop','orig_ctop2',...
        'sarta_clr_co2_used','sarta_xclr_co2_used','sarta_cld_co2_used','sarta_rclearcalc'};
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end
	  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% these are from PCRTM
clear list
list = {'rad_allsky_std','rcalc_std','pcrtm_lvlODice','pcrtm_lvlODwater','pcrtm_iceODX','pcrtm_waterODX',...
        'pcrtm_iceOD','pcrtm_iceDME','pcrtm_iceCTOP','pcrtm_waterOD','pcrtm_waterDME','pcrtm_waterCTOP',...
	'pcrtm_co2_used','ncol','overlap','rad_allsky','rad_clrsky'};

for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end
fprintf(1,'blew away %2i fields \n',iF);
