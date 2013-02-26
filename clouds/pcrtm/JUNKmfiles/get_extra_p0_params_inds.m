p0ALL.robs1(:,inds)       = profRX.robs1;
p0ALL.sarta_clear(:,inds) = profRX.rcalc;
p0ALL.rad_allsky(:,inds)  = rad_allsky(h.ichan,:)*1000;
p0ALL.rad_clrsky(:,inds)  = rad_clrsky(h.ichan,:)*1000;

p0ALL.ncol(inds)                = ones(size(p0ALL.cpsize(inds)))*ncol0;
p0ALL.sarta_clr_co2_used(inds)  = ppmvLAY(40,:);
p0ALL.pcrtm_co2_used(inds)      = co2;
p0ALL.overlap(inds)             = ones(size(p0ALL.cpsize(inds)))*overlap;

p0ALL.rlat(inds)         = p.rlat;
p0ALL.rlon(inds)         = p.rlon;
p0ALL.landfrac(inds)     = p.landfrac;
p0ALL.solzen(inds)       = p.solzen;
p0ALL.stemp(inds)        = p.stemp;

p0ALL.ctype(inds)        = p.ctype;
p0ALL.cfrac(inds)        = p.cfrac;
p0ALL.cngwat(inds)       = p.cngwat;
p0ALL.cpsize(inds)       = p.cpsize;
p0ALL.cprtop(inds)       = p.cprtop;

p0ALL.ctype2(inds)       = p.ctype2;
p0ALL.cfrac2(inds)       = p.cfrac2;
p0ALL.cngwat2(inds)      = p.cngwat2;
p0ALL.cpsize2(inds)      = p.cpsize2;
p0ALL.cprtop2(inds)      = p.cprtop2;

p0ALL.pcrtm_iceOD(inds)     = tmpjunk.totalODice;
p0ALL.pcrtm_iceDME(inds)    = tmpjunk.meanDMEice;
p0ALL.pcrtm_iceCTOP(inds)   = tmpjunk.maxCTOPice;
p0ALL.pcrtm_waterOD(inds)   = tmpjunk.totalODwater;
p0ALL.pcrtm_waterDME(inds)  = tmpjunk.meanDMEwater;
p0ALL.pcrtm_waterCTOP(inds) = tmpjunk.maxCTOPwater;
