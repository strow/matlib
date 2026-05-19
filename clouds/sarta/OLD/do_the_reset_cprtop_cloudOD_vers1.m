function p1 = do_the_reset_cprtop_cloudOD_vers1(p0,p1,pICE,pWATER);

ODice   = p1.sarta_lvlODice;
ODwater = p1.sarta_lvlODwater;

iDebug = -1;
if iDebug > 0
  junky = 105;
  p10 = p1;
  p00 = p0;
  disp('blah.m')
  [p1.cprtop(junky) p1.cprbot(junky) p1.cprtop2(junky) p1.cprbot2(junky)]
end

ice = find(p1.ctype == 201 & p1.ctype2 ~= 201 & pICE > 0 & p1.sarta_iceOD_warn < 0);
for oo = 1 : length(ice)
  od = ODice(:,ice(oo)); 
  plevs = p1.plevs(:,ice(oo));
  peakindex = p1.sarta_index_wgtpeakI(ice(oo));  %% use one index above this for actual top
  if peakindex < 2
    peakindex = 2;
  end
  cldthick = p1.cprbot(ice(oo)) - p1.cprtop(ice(oo));
  p1.cprtop(ice(oo)) = max(10,p1.plevs(peakindex-1));
  p1.cprbot(ice(oo)) = p1.cprtop(ice(oo)) + cldthick;
end
if iDebug > 0
  [p1.cprtop(junky) p1.cprbot(junky) p1.cprtop2(junky) p1.cprbot2(junky)]
end

ice = find(p1.ctype2 == 201 & p1.ctype ~= 201 & pICE > 0 & p1.sarta_iceOD_warn < 0);
for oo = 1 : length(ice)
  od = ODice(:,ice(oo)); 
  plevs = p1.plevs(:,ice(oo));
  peakindex = p1.sarta_index_wgtpeakI(ice(oo));  %% use one index above this for actual top
  if peakindex < 2
    peakindex = 2;
  end
  cldthick = p1.cprbot2(ice(oo)) - p1.cprtop2(ice(oo));
  p1.cprtop2(ice(oo)) = max(10,p1.plevs(peakindex-1));
  p1.cprbot2(ice(oo)) = p1.cprtop2(ice(oo)) + cldthick;
end
if iDebug > 0
  [p1.cprtop(junky) p1.cprbot(junky) p1.cprtop2(junky) p1.cprbot2(junky)]
end

ice = find(p1.ctype2 == 201 & p1.ctype == 201 & pICE > 0 & p1.sarta_iceOD_warn < 0);
for oo = 1 : length(ice)
  od = ODice(:,ice(oo)); 
  plevs = p1.plevs(:,ice(oo));
  peakindex = p1.sarta_index_wgtpeakI(ice(oo));  %% use one index above this for actual top
  if peakindex < 2
    peakindex = 2;
  end
  cldthick = p1.cprbot(ice(oo)) - p1.cprtop(ice(oo));
  p1.cprtop(ice(oo)) = max(10,p1.plevs(peakindex-1));
  p1.cprbot(ice(oo)) = p1.cprtop(ice(oo)) + cldthick;
end
if iDebug > 0
  [p1.cprtop(junky) p1.cprbot(junky) p1.cprtop2(junky) p1.cprbot2(junky)]
end

%% this is water

water = find(p1.ctype == 101 & p1.ctype2 ~= 101 & pWATER > 0 & p1.sarta_waterOD_warn < 0);
for oo = 1 : length(water)
  od = ODwater(:,water(oo)); 
  plevs = p1.plevs(:,water(oo));
  peakindex = p1.sarta_index_wgtpeakW(water(oo));  %% use one index above this for actual top
  if peakindex < 2
    peakindex = 2;
  end
  cldthick = p1.cprbot(water(oo)) - p1.cprtop(water(oo));
  p1.cprtop(water(oo)) = max(10,p1.plevs(peakindex-1));
  p1.cprbot(water(oo)) = p1.cprtop(water(oo)) + cldthick;
end
if iDebug > 0
  [p1.cprtop(junky) p1.cprbot(junky) p1.cprtop2(junky) p1.cprbot2(junky)]
end

water = find(p1.ctype2 == 101 & p1.ctype ~= 101 & pWATER > 0 & p1.sarta_waterOD_warn < 0);
for oo = 1 : length(water)
  od = ODwater(:,water(oo)); 
  plevs = p1.plevs(:,water(oo));
  peakindex = p1.sarta_index_wgtpeakW(water(oo));  %% use one index above this for actual top
  if peakindex < 2
    peakindex = 2;
  end
  cldthick = p1.cprbot(water(oo)) - p1.cprtop(water(oo));
  p1.cprtop2(water(oo)) = max(10,p1.plevs(peakindex-1));
  p1.cprbot2(water(oo)) = p1.cprtop2(water(oo)) + cldthick;
end
if iDebug > 0
  [p1.cprtop(junky) p1.cprbot(junky) p1.cprtop2(junky) p1.cprbot2(junky)]
end

water = find(p1.ctype2 == 101 & p1.ctype == 101 & pWATER > 0 & p1.sarta_waterOD_warn < 0);
for oo = 1 : length(water)
  od = ODwater(:,water(oo)); 
  plevs = p1.plevs(:,water(oo));
  peakindex = p1.sarta_index_wgtpeakW(water(oo));  %% use one index above this for actual top
  if peakindex < 2
    peakindex = 2;
  end
  cldthick = p1.cprbot(water(oo)) - p1.cprtop(water(oo));
  p1.cprtop(water(oo)) = max(10,p1.plevs(peakindex-1));
  p1.cprbot(water(oo)) = p1.cprtop(water(oo)) + cldthick;
end

if iDebug > 0
  [p1.cprtop(junky) p1.cprbot(junky) p1.cprtop2(junky) p1.cprbot2(junky)]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now do the corrections based on cc/ciwc and cc/clwc
ice1 = find(p1.ctype2 == 201 & p1.ctype == 201 & p1.sarta_iceOD_warn == 1);
  p1.cngwat(ice1) = 0;  
  p1.cngwat2(ice1) = 0; 
  p1.cfrac(ice1) = 0;  
  p1.cfrac2(ice1) = 0; 
  p1.cfrac12(ice1) = 0; 
ice2 = find(p1.ctype2 == 201 & p1.ctype ~= 201& p1.sarta_iceOD_warn == 1);
  p1.cngwat2(ice2) = 0; 
  p1.cfrac2(ice2) = 0; 
  p1.cfrac12(ice2) = 0; 
ice3 = find(p1.ctype == 201 & p1.ctype2 ~= 201 & p1.sarta_iceOD_warn == 1);
  p1.cngwat(ice3) = 0;  
  p1.cfrac(ice3) = 0;  
  p1.cfrac12(ice3) = 0; 
icebad = length([ice1 ice2 ice3]); 

water1 = find(p1.ctype2 == 101 & p1.ctype == 101 & p1.sarta_waterOD_warn == 1);
  p1.cngwat(water1) = 0;  
  p1.cngwat2(water1) = 0; 
  p1.cfrac(water1) = 0;  
  p1.cfrac2(water1) = 0; 
  p1.cfrac12(water1) = 0; 
water2 = find(p1.ctype2 == 101 & p1.ctype ~= 101& p1.sarta_waterOD_warn == 1);
  p1.cngwat2(water2) = 0; 
  p1.cfrac2(water2) = 0; 
  p1.cfrac12(water2) = 0; 
water3 = find(p1.ctype == 101 & p1.ctype2 ~= 101 & p1.sarta_waterOD_warn == 1);
  p1.cngwat(water3) = 0;  
  p1.cfrac(water3) = 0;  
  p1.cfrac12(water3) = 0; 
waterbad = length([water1 water2 water3]);

fprintf(1,'  looking at coincidence of ciwc/cc and clwc/cc, reset %5i bad ice and %5i bad water profs \n',icebad,waterbad)