%if isfield(p0ALL,'robs1')
%  p1ALL.robs1(:,inds)       = p0ALL.robs1(:,inds);
%end

if run_sarta.clear > 0
  p1ALL.sarta_clear(:,inds)       = profRX.rcalc;
  p1ALL.sarta_clr_co2_used(inds)  = ppmvLAY(40,:);
end

if run_sarta.cloud > 0
  p1ALL.sarta_cloud(:,inds)       = profRX2.rcalc;
  p1ALL.sarta_cld_co2_used(inds)  = ppmvLAY2(40,:);

  if isfield(p2junk,'sarta_lvlZ')
    p1ALL.sarta_lvlZ(:,inds) = p2junk.sarta_lvlZ;
  end

  if isfield(p2junk,'icecldX')
    p1ALL.icecldX(:,inds) = p2junk.icecldX;
  end
  if isfield(p2junk,'icecldY')
    p1ALL.icecldY(:,inds) = p2junk.icecldY;
  end
  if isfield(p2junk,'sarta_lvlDMEice')
    p1ALL.sarta_lvlDMEice(:,inds) = p2junk.sarta_lvlDMEice;
  end
  if isfield(p2junk,'sarta_lvlODice')
    p1ALL.sarta_lvlODice(:,inds) = p2junk.sarta_lvlODice;
  end
  if isfield(p2junk,'sarta_lvl_iceOD_1')
    p1ALL.sarta_lvl_iceOD_1(inds) = p2junk.sarta_lvl_iceOD_1;
  end
  if isfield(p2junk,'sarta_wgtpeakI')
    p1ALL.sarta_wgtpeakI(inds) = p2junk.sarta_wgtpeakI;
  end

  if isfield(p2junk,'watercldX')
    p1ALL.watercldX(:,inds) = p2junk.watercldX;
  end
  if isfield(p2junk,'watercldY')
    p1ALL.watercldY(:,inds) = p2junk.watercldY;
  end
  if isfield(p2junk,'sarta_lvlDMEwater')
    p1ALL.sarta_lvlDMEwater(:,inds) = p2junk.sarta_lvlDMEwater;
  end
  if isfield(p2junk,'sarta_lvlODwater')
    p1ALL.sarta_lvlODwater(:,inds) = p2junk.sarta_lvlODwater;
  end
  if isfield(p2junk,'sarta_lvl_waterOD_1')
    p1ALL.sarta_lvl_waterOD_1(inds) = p2junk.sarta_lvl_waterOD_1;
  end
  if isfield(p2junk,'sarta_wgtpeakI')
    p1ALL.sarta_wgtpeakW(inds) = p2junk.sarta_wgtpeakW;
  end

end

p1ALL.rad_allsky(:,inds)  = rad_allsky(h.ichan,:)*1000;
p1ALL.rad_allsky_std(:,inds)  = rad_allsky_std(h.ichan,:)*1000;

p1ALL.rad_clrsky(:,inds)  = rad_clrsky(h.ichan,:)*1000;

p1ALL.ncol(inds)                = ones(size(p0ALL.stemp(inds)))*ncol0;
p1ALL.pcrtm_co2_used(inds)      = co2;
p1ALL.overlap(inds)             = ones(size(p0ALL.stemp(inds)))*overlap;

p1ALL.rlat(inds)         = p.rlat;
p1ALL.rlon(inds)         = p.rlon;
p1ALL.landfrac(inds)     = p.landfrac;
p1ALL.solzen(inds)       = p.solzen;
p1ALL.stemp(inds)        = p.stemp;

%% these are the SARTA two slab CLOUD parameters
if isfield(p,'ctype')
  p1ALL.ctype(inds)        = p.ctype;
end
if isfield(p,'cfrac')
  p1ALL.cfrac(inds)        = p.cfrac;
end
if isfield(p,'cngwat')
  p1ALL.cngwat(inds)       = p.cngwat;
end
if isfield(p,'cpsize')
  p1ALL.cpsize(inds)       = p.cpsize;
end
if isfield(p,'cprtop')
  p1ALL.cprtop(inds)       = p.cprtop;
end
if isfield(p,'cprbot')
  p1ALL.cprbot(inds)       = p.cprbot;
end

if isfield(p,'ctype2')
  p1ALL.ctype2(inds)        = p.ctype2;
end
if isfield(p,'cfrac2')
  p1ALL.cfrac2(inds)        = p.cfrac2;
end
if isfield(p,'cngwat2')
  p1ALL.cngwat2(inds)       = p.cngwat2;
end
if isfield(p,'cpsize2')
  p1ALL.cpsize2(inds)       = p.cpsize2;
end
if isfield(p,'cprtop2')
  p1ALL.cprtop2(inds)       = p.cprtop2;
end
if isfield(p,'cprbot2')
  p1ALL.cprbot2(inds)       = p.cprbot2;
end

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

