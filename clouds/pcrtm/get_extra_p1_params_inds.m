p1ALL.robs1(:,inds)       = profRX.robs1;

if run_sarta.clear > 0
  p1ALL.sarta_clear(:,inds) = profRX.rcalc;
  p1ALL.sarta_clr_co2_used(inds)  = ppmvLAY(40,:);
end

if run_sarta.cloud > 0
  p1ALL.sarta_cloud(:,inds) = profRX2.rcalc;
  p1ALL.sarta_cld_co2_used(inds)  = ppmvLAY2(40,:);
end

p1ALL.rad_allsky(:,inds)  = rad_allsky(h.ichan,:)*1000;
p1ALL.rad_clrsky(:,inds)  = rad_clrsky(h.ichan,:)*1000;

p1ALL.ncol(inds)                = ones(size(p0ALL.cpsize(inds)))*ncol0;
p1ALL.pcrtm_co2_used(inds)      = co2;
p1ALL.overlap(inds)             = ones(size(p0ALL.cpsize(inds)))*overlap;

p1ALL.rlat(inds)         = p.rlat;
p1ALL.rlon(inds)         = p.rlon;
p1ALL.landfrac(inds)     = p.landfrac;
p1ALL.solzen(inds)       = p.solzen;
p1ALL.stemp(inds)        = p.stemp;

p1ALL.ctype(inds)        = p.ctype;
p1ALL.cfrac(inds)        = p.cfrac;
p1ALL.cngwat(inds)       = p.cngwat;
p1ALL.cpsize(inds)       = p.cpsize;
p1ALL.cprtop(inds)       = p.cprtop;

p1ALL.ctype2(inds)       = p.ctype2;
p1ALL.cfrac2(inds)       = p.cfrac2;
p1ALL.cngwat2(inds)      = p.cngwat2;
p1ALL.cpsize2(inds)      = p.cpsize2;
p1ALL.cprtop2(inds)      = p.cprtop2;

p1ALL.pcrtm_iceOD(inds)     = tmpjunk.totalODice;
p1ALL.pcrtm_iceDME(inds)    = tmpjunk.meanDMEice;
p1ALL.pcrtm_iceCTOP(inds)   = tmpjunk.maxCTOPice;
p1ALL.pcrtm_waterOD(inds)   = tmpjunk.totalODwater;
p1ALL.pcrtm_waterDME(inds)  = tmpjunk.meanDMEwater;
p1ALL.pcrtm_waterCTOP(inds) = tmpjunk.maxCTOPwater;
