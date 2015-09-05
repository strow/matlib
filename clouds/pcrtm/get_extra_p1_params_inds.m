%% update PCRTM output

p1ALL.rad_allsky(:,inds)     = rad_allsky(h.ichan,:)*1000;
p1ALL.rad_allsky_std(:,inds) = rad_allsky_std(h.ichan,:)*1000;
p1ALL.rad_clrsky(:,inds)     = rad_clrsky(h.ichan,:)*1000;

p1ALL.ncol(inds)           = ones(size(p0ALL.stemp(inds)))*ncol0;
p1ALL.pcrtm_co2_used(inds) = co2;
p1ALL.overlap(inds)        = ones(size(p0ALL.stemp(inds)))*overlap;

p1ALL.rlat(inds)         = p.rlat;
p1ALL.rlon(inds)         = p.rlon;
p1ALL.landfrac(inds)     = p.landfrac;
p1ALL.solzen(inds)       = p.solzen;
p1ALL.stemp(inds)        = p.stemp;

%% these are the PCRTM MEAN CLOUD column parameters, which is basically : sub_opt = ice(high) + water(low)
p1ALL.pcrtm_iceOD(inds)     = tmpjunk.totalODice;
p1ALL.pcrtm_iceDME(inds)    = tmpjunk.meanDMEice;
p1ALL.pcrtm_iceCTOP(inds)   = tmpjunk.maxCTOPice;
p1ALL.pcrtm_waterOD(inds)   = tmpjunk.totalODwater;
p1ALL.pcrtm_waterDME(inds)  = tmpjunk.meanDMEwater;
p1ALL.pcrtm_waterCTOP(inds) = tmpjunk.maxCTOPwater;

%% these are the PCRTM CLOUD profile parameters, and should compare to each other
p1ALL.pcrtm_lvlODice(:,inds)   = tmpjunk.lvlODice;
p1ALL.pcrtm_lvlODwater(:,inds) = tmpjunk.lvlODwater;
p1ALL.pcrtm_iceODX(inds)       = tmpjunk.totalODiceX;
p1ALL.pcrtm_waterODX(inds)     = tmpjunk.totalODwaterX;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% update SARTA clear output
if run_sarta.clear > 0
  p1ALL.sarta_rclearcalc(:,inds)  = profRX.rcalc;  %% this is with PCRTM co2 profile
  p1ALL.sarta_clear(:,inds)       = profRX.rcalc;  %% this is with PCRTM co2 profile
  p1ALL.sarta_clr_co2_used(inds)  = ppmvLAY(40,:);
  if exist('xprofRX')
    p1ALL.sarta_xclear(:,inds)       = xprofRX.rcalc;  %% this is with SARTA co2 profile
    p1ALL.sarta_xclr_co2_used(inds)  = xppmvLAY(40,:);     
    disp('   >>> added both clear calcs (default SARTA profile and PCRTM profile) ...')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% update SARTA cloud output

if run_sarta.cloud > 0
  p1ALL.sarta_rclearcalc(:,inds)  = profRX.rcalc;   
  p1ALL.sarta_cloud(:,inds)       = profRX2.rcalc;
  p1ALL.sarta_cld_co2_used(inds)  = ppmvLAY2(40,:);

  if isfield(p2junk,'sarta_lvlZ')
    p1ALL.sarta_lvlZ(:,inds) = p2junk.sarta_lvlZ;
  end

  listjunk = {'ctype','cfrac','cngwat','cpsize','cprtop','cprbot','ctype2','cfrac2','cngwat2','cpsize2','cprtop2','cprbot2','cfrac12'};
  for iijunk = 1 : length(listjunk)
    fieldjunk = listjunk{iijunk};
    if isfield(p2junk,fieldjunk)
      strjunk = ['p1ALL.' fieldjunk '(inds) = p2junk.' fieldjunk ';'];
      %fprintf(1,'%s \n',strjunk)
      eval(strjunk)
    end
  end
  clear listjunk iijunk fieldjunk
  
  listjunk = {'icecldX','icecldY','watercldX','watercldY','sarta_lvl_iceOD_1','sarta_wgtpeakI','sarta_lvl_waterOD_1','sarta_wgtpeakW'};
  for iijunk = 1 : length(listjunk)
    fieldjunk = listjunk{iijunk};
    if isfield(p2junk,fieldjunk)
      strjunk = ['p1ALL.' fieldjunk '(inds) = p2junk.' fieldjunk ';'];
      %fprintf(1,'%s \n',strjunk)
      eval(strjunk)
    end
  end
  clear listjunk iijunk fieldjunk
  
  listjunk = {'sarta_lvlDMEice','sarta_lvlODice','sarta_lvlDMEwater','sarta_lvlODwater'};
  for iijunk = 1 : length(listjunk)
    fieldjunk = listjunk{iijunk};
    if isfield(p2junk,fieldjunk)
      strjunk = ['p1ALL.' fieldjunk '(:,inds) = p2junk.' fieldjunk ';'];
      %fprintf(1,'%s \n',strjunk)      
      eval(strjunk)
    end
  end
  clear listjunk iijunk fieldjunk  
  
elseif run_sarta.cloud < 0

  %% these are the orig input SARTA two slab CLOUD parameters, if they exist
  %% keep as is
  listjunk = {'ctype','cfrac','cngwat','cpsize','cprtop','cprbot','ctype2','cfrac2','cngwat2','cpsize2','cprtop2','cprbot2','cfrac12'};
  for iijunk = 1 : length(listjunk)
    fieldjunk = listjunk{iijunk};
    if isfield(p,fieldjunk)
      strjunk = ['p1ALL.' fieldjunk '(inds) = p.' fieldjunk ';'];
      eval(str)
    end
  end
  clear listjunk iijunk fieldjunk

end