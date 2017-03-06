function p1 = reset_cprtop(p0);

p1 = p0;
p1.orig_ctop  = p1.cprtop;
p1.orig_ctop2 = p1.cprtop2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ecmwfcld2sartacld.m -- > new_style_smooth_cc_ciwc_clwc_to_water_ice_profile --> cloud_mean_press --> set icecldY,watercldY
%% from cloud_mean_press we get
%%   xcumsum = run_sarta.cumsum so can be -9999,-1,0-1,1-9998,9999
%%   aa.icecldX,aa.watercldX = mean(CIWC), mean(CLWC) pressure level
%%   aa.icecldY,aa.watercldY = pressure level where normalized CIWC/CLWC exceed xcumsum if 0 < xcumsum < 1
%%                             else set to 1200 mb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is ice

ice = find(p1.ctype == 201 & p1.ctype2 ~= 201);
pcenter_ice(ice) = (p0.cprtop(ice) + p0.cprbot(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > p1.icecldY(ice)); 
if length(oo) > 0
  delta = pcenter_ice(ice) - p1.icecldY(ice);
  p1.cprtop(ice(oo)) = p1.cprtop(ice(oo)) - delta(oo);
  p1.cprbot(ice(oo)) = p1.cprbot(ice(oo)) - delta(oo);
end

ice = find(p1.ctype2 == 201 & p1.ctype ~= 201);
pcenter_ice(ice) = (p0.cprtop2(ice) + p0.cprbot2(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > p1.icecldY(ice)); 
if length(oo) > 0
  delta = pcenter_ice(ice) - p1.icecldY(ice);
  p1.cprtop2(ice(oo)) = p1.cprtop2(ice(oo)) - delta(oo);
  p1.cprbot2(ice(oo)) = p1.cprbot2(ice(oo)) - delta(oo);
end

ice = find(p1.ctype2 == 201 & p1.ctype == 201);
pcenter_ice(ice) = (p0.cprtop(ice) + p0.cprbot(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > p1.icecldY(ice)); 
if length(oo) > 0
  delta = pcenter_ice(ice) - p1.icecldY(ice);
  p1.cprtop(ice(oo)) = p1.cprtop(ice(oo)) - delta(oo);
  p1.cprbot(ice(oo)) = p1.cprbot(ice(oo)) - delta(oo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is water

water = find(p1.ctype == 101 & p1.ctype2 ~= 101);
pcenter_water(water) = (p0.cprtop(water) + p0.cprbot(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > p1.watercldY(water)); 
if length(oo) > 0
  delta = pcenter_water(water) - p1.watercldY(water);
  p1.cprtop(water(oo)) = p1.cprtop(water(oo)) - delta(oo);
  p1.cprbot(water(oo)) = p1.cprbot(water(oo)) - delta(oo);
end

water = find(p1.ctype2 == 101 & p1.ctype ~= 101);
pcenter_water(water) = (p0.cprtop2(water) + p0.cprbot2(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > p1.watercldY(water)); 
if length(oo) > 0
  delta = pcenter_water(water) - p1.watercldY(water);
  p1.cprtop2(water(oo)) = p1.cprtop2(water(oo)) - delta(oo);
  p1.cprbot2(water(oo)) = p1.cprbot2(water(oo)) - delta(oo);
end

water = find(p1.ctype2 == 101 & p1.ctype == 101);
pcenter_water(water) = (p0.cprtop(water) + p0.cprbot(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > p1.watercldY(water)); 
if length(oo) > 0
  delta = pcenter_water(water) - p1.watercldY(water);
  p1.cprtop(water(oo)) = p1.cprtop(water(oo)) - delta(oo);
  p1.cprbot(water(oo)) = p1.cprbot(water(oo)) - delta(oo);
end

iDebug = -1;
if iDebug > 0
  disp('resetting cloud top/bot in reset_cprtop.m')
  com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false)   
  keyboard
  plot(1:length(p0.stemp),p0.cprtop,'b',1:length(p0.stemp),p0.cprbot,'c',1:length(p0.stemp),p1.cprtop,'r',1:length(p0.stemp),p1.cprbot,'m')
%  plot(1:length(p0.stemp),(p0.cprtop-p0.cprbot) - (p1.cprtop-p1.cprbot),'o-')
end
