function p1 = do_the_reset_cprtop_cloudOD_vers0(p0,p1,pICE,pWATER);

%% this is ice

ice = find(p1.ctype == 201 & p1.ctype2 ~= 201);
pcenter_ice(ice) = (p0.cprtop(ice) + p0.cprbot(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > pICE(ice) & pICE(ice) > 0); 
if length(oo) > 0
  delta = pcenter_ice(ice) - pICE(ice);
  p1.cprtop(ice(oo)) = p1.cprtop(ice(oo)) - delta(oo);
  p1.cprbot(ice(oo)) = p1.cprbot(ice(oo)) - delta(oo);
end

ice = find(p1.ctype2 == 201 & p1.ctype ~= 201);
pcenter_ice(ice) = (p0.cprtop2(ice) + p0.cprbot2(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > pICE(ice) & pICE(ice) > 0); 
if length(oo) > 0
  delta = pcenter_ice(ice) - pICE(ice);
  p1.cprtop2(ice(oo)) = p1.cprtop2(ice(oo)) - delta(oo);
  p1.cprbot2(ice(oo)) = p1.cprbot2(ice(oo)) - delta(oo);
end

ice = find(p1.ctype2 == 201 & p1.ctype == 201);
pcenter_ice(ice) = (p0.cprtop(ice) + p0.cprbot(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > pICE(ice) & pICE(ice) > 0);  
if length(oo) > 0
  delta = pcenter_ice(ice) - pICE(ice);
  p1.cprtop(ice(oo)) = p1.cprtop(ice(oo)) - delta(oo);
  p1.cprbot(ice(oo)) = p1.cprbot(ice(oo)) - delta(oo);
end

%% this is water

water = find(p1.ctype == 101 & p1.ctype2 ~= 101);
pcenter_water(water) = (p0.cprtop(water) + p0.cprbot(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > pWATER(water) & pWATER(water) > 0);
if length(oo) > 0
  delta = pcenter_water(water) - pWATER(water);
  p1.cprtop(water(oo)) = p1.cprtop(water(oo)) - delta(oo);
  p1.cprbot(water(oo)) = p1.cprbot(water(oo)) - delta(oo);
end

water = find(p1.ctype2 == 101 & p1.ctype ~= 101);
pcenter_water(water) = (p0.cprtop2(water) + p0.cprbot2(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > pWATER(water) & pWATER(water) > 0);
if length(oo) > 0
  delta = pcenter_water(water) - pWATER(water);
  p1.cprtop2(water(oo)) = p1.cprtop2(water(oo)) - delta(oo);
  p1.cprbot2(water(oo)) = p1.cprbot2(water(oo)) - delta(oo);
end

water = find(p1.ctype2 == 101 & p1.ctype == 101);
pcenter_water(water) = (p0.cprtop(water) + p0.cprbot(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > pWATER(water) & pWATER(water) > 0);
if length(oo) > 0
  delta = pcenter_water(water) - pWATER(water);
  p1.cprtop(water(oo)) = p1.cprtop(water(oo)) - delta(oo);
  p1.cprbot(water(oo)) = p1.cprbot(water(oo)) - delta(oo);
end
