p0.robs1       = profRX.robs1;
p0.sarta_clear = profRX.rcalc;
p0.rad_allsky  = rad_allsky(h.ichan,:)*1000;
p0.rad_clrsky  = rad_clrsky(h.ichan,:)*1000;
p0.ncol        = ones(size(p0.cpsize))*ncol0;
p0.sarta_clr_co2_used  = ppmvLAY(40,:);
p0.pcrtm_co2_used      = co2;
p0.overlap             = ones(size(p0.cpsize))*overlap;

p0.rlat         = p.rlat;
p0.rlon         = p.rlon;
p0.landfrac     = p.landfrac;
p0.solzen       = p.solzen;
p0.stemp        = p.stemp;

p0.ctype        = p.ctype;
p0.cfrac        = p.cfrac;
p0.cngwat       = p.cngwat;
p0.cpsize       = p.cpsize;
p0.cprtop       = p.cprtop;

p0.ctype2       = p.ctype2;
p0.cfrac2       = p.cfrac2;
p0.cngwat2      = p.cngwat2;
p0.cpsize2      = p.cpsize2;
p0.cprtop2      = p.cprtop2;

p0.pcrtm_iceOD   = tmpjunk.totalODice;
p0.pcrtm_iceDME  = tmpjunk.meanDMEice;
p0.pcrtm_iceCTOP = tmpjunk.maxCTOPice;
p0.pcrtm_waterOD   = tmpjunk.totalODwater;
p0.pcrtm_waterDME  = tmpjunk.meanDMEwater;
p0.pcrtm_waterCTOP = tmpjunk.maxCTOPwater;
