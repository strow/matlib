function p1 = reset_cprtop_cloudOD(p0);

%% computes cloud ODs based on formulas given by Xianglei and Xiuhong
%% see PCRTM_compute_for_AIRS_spectra.m

p1 = p0;
p1.orig_ctop  = p1.cprtop;
p1.orig_ctop2 = p1.cprtop2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 287.05;
g = 9.806;

% coefficents from S-C Ou, K-N. Liou, Atmospheric Research
% 35(1995):127-138.
% for computing ice cloud effective size
c0 = 326.3;
c1 = 12.42;
c2 = 0.197;
c3 = 0.0012;

for ii = 1 : length(p0.stemp)
  ciwc = p0.ciwc(:,ii);
  clwc = p0.clwc(:,ii);
  cc   = p0.cc(:,ii);
  ptemp = p0.ptemp(:,ii);
  gas_1 = p0.gas_1(:,ii);
  press = p0.plevs(:,ii);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %% code works assming press(1) < press(2) etc ie TOA is at LVL1, GND = LVLN
  if press(1) > press(2)
    [Y,I] = sort(press);
    clwc  = clwc(I);
    ciwc  = ciwc(I);
    cc    = cc(I);
    ptemp = ptemp(I);
    gas_1 = gas_1(I);
    press = press(I);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  tcld = ptemp - 273.16;                   %% change to deg C
  boo = find(tcld < -50); tcld(boo) = -50; %% set minimum as -50 C
  boo = find(tcld > -25); tcld(boo) = -25; %% set maximum as -50 C

  cldde_ice = c0 + c1 * tcld + c2 * tcld.^2 + c3 * tcld.^3;

  tav = (ptemp(1:end-1)+ptemp(2:end))/2; tav(end) = ptemp(end);
  pav = (press(1:end-1)+press(2:end))/2; pav(end) = press(end);
  scaleH = R*tav/g/1000;    %% in km
  dz = scaleH .* log(press(2:end)./press(1:end-1)); dz(length(dz)+1) = 0;
  Z = cumsum(dz);
  diffZ = abs(diff(Z)); diffZ(length(diffZ)+1) = diffZ(length(diffZ));

  % compute ice cloud optical depth from Ebert and Curry (1992, J. Geophys. Res.,  
  qi = ciwc ./ ptemp .* press *100/R * 1e3;  %%change IWC from kg/kg to g/m3  
  iceOD = (0.003448 + 2.431./cldde_ice) .*qi ./ cc .* diffZ *1e3;

  bad = find(isnan(cldde_ice)); cldde_ice(bad) = -9999;
  p1.sarta_lvlDMEice(:,ii) = cldde_ice;
  bad = find(isnan(iceOD)); iceOD(bad) = 0;
  p1.sarta_lvlODice(:,ii) = iceOD;
  p1.sarta_lvlZ(:,ii) = Z;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  cldde_liq = 20 * ones(size(press));

  qw = clwc ./ ptemp .* press *100/R *1e3;  %change liquid water content from kg/kg to g/m^3
  % ECMWF technical report, Klein S. A., 1999
  waterOD = 3 * qw ./ cldde_liq ./cc .*diffZ  *1e3;

  bad = find(isnan(cldde_liq)); cldde_liq(bad) = -9999;
  p1.sarta_lvlDMEwater(:,ii) = cldde_liq;
  bad = find(isnan(waterOD)); waterOD(bad) = 0;
  p1.sarta_lvlODwater(:,ii) = waterOD;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  icesum = cumsum(iceOD);
  watersum = cumsum(waterOD);

  od_ice1_top = find(icesum > 1,1);
  if length(od_ice1_top) > 0
    p1.sarta_lvl_iceOD_1(ii) = press(od_ice1_top);
  else
    oo = find(iceOD > 0);
    if length(oo) > 0
      %% cumulative OD < 1, put as low as possible
      oo = oo(end);
      od_ice1_top = oo;
      p1.sarta_lvl_iceOD_1(ii) = press(oo);
    else
      od_ice1_top = NaN;
      p1.sarta_lvl_iceOD_1(ii) = -9999;
    end
  end

  od_water1_top = find(watersum > 1,1);
  if length(od_water1_top) > 0
    p1.sarta_lvl_waterOD_1(ii) = press(od_water1_top);
  else
    oo = find(waterOD > 0);
    if length(oo) > 0
      %% cumulative OD < 1, put as low as possible
      oo = oo(end);
      od_water1_top = oo;
      p1.sarta_lvl_waterOD_1(ii) = press(oo);
    else
      od_water1_top = NaN;
      p1.sarta_lvl_waterOD_1(ii) = -9999;
    end
  end

%{
  plot(1:length(press),iceOD/nansum(iceOD),'b',1:length(press),waterOD/nansum(waterOD),'r',...
       od_ice1_top,0.1,'bo',       od_water1_top,0.1,'ro')
  plot(1:length(press),iceOD,'b',1:length(press),waterOD,'r',...
       od_ice1_top,0.1,'bx',       od_water1_top,0.2,'ro')
  title(num2str(ii));
  pause(0.1)
%}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is basically "reset_cprtop"
%% except we replaced 
%%   icecldY   with sarta_lvl_iceOD_1
%%   watercldY with sarta_lvl_waterOD_1

%% this is ice

ice = find(p1.ctype == 201 & p1.ctype2 ~= 201);
pcenter_ice(ice) = (p0.cprtop(ice) + p0.cprbot(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > p1.sarta_lvl_iceOD_1(ice) & p1.sarta_lvl_iceOD_1(ice) > 0); 
if length(oo) > 0
  delta = pcenter_ice(ice) - p1.sarta_lvl_iceOD_1(ice);
  p1.cprtop(ice(oo)) = p1.cprtop(ice(oo)) - delta(oo);
  p1.cprbot(ice(oo)) = p1.cprbot(ice(oo)) - delta(oo);
end

ice = find(p1.ctype2 == 201 & p1.ctype ~= 201);
pcenter_ice(ice) = (p0.cprtop2(ice) + p0.cprbot2(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > p1.sarta_lvl_iceOD_1(ice) & p1.sarta_lvl_iceOD_1(ice) > 0); 
if length(oo) > 0
  delta = pcenter_ice(ice) - p1.sarta_lvl_iceOD_1(ice);
  p1.cprtop2(ice(oo)) = p1.cprtop2(ice(oo)) - delta(oo);
  p1.cprbot2(ice(oo)) = p1.cprbot2(ice(oo)) - delta(oo);
end

ice = find(p1.ctype2 == 201 & p1.ctype == 201);
pcenter_ice(ice) = (p0.cprtop(ice) + p0.cprbot(ice))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_ice(ice) > p1.sarta_lvl_iceOD_1(ice) & p1.sarta_lvl_iceOD_1(ice) > 0);  
if length(oo) > 0
  delta = pcenter_ice(ice) - p1.sarta_lvl_iceOD_1(ice);
  p1.cprtop(ice(oo)) = p1.cprtop(ice(oo)) - delta(oo);
  p1.cprbot(ice(oo)) = p1.cprbot(ice(oo)) - delta(oo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is water

water = find(p1.ctype == 101 & p1.ctype2 ~= 101);
pcenter_water(water) = (p0.cprtop(water) + p0.cprbot(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > p1.sarta_lvl_waterOD_1(water) & p1.sarta_lvl_waterOD_1(water) > 0);
if length(oo) > 0
  delta = pcenter_water(water) - p1.sarta_lvl_waterOD_1(water);
  p1.cprtop(water(oo)) = p1.cprtop(water(oo)) - delta(oo);
  p1.cprbot(water(oo)) = p1.cprbot(water(oo)) - delta(oo);
end

water = find(p1.ctype2 == 101 & p1.ctype ~= 101);
pcenter_water(water) = (p0.cprtop2(water) + p0.cprbot2(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > p1.sarta_lvl_waterOD_1(water) & p1.sarta_lvl_waterOD_1(water) > 0);
if length(oo) > 0
  delta = pcenter_water(water) - p1.sarta_lvl_waterOD_1(water);
  p1.cprtop2(water(oo)) = p1.cprtop2(water(oo)) - delta(oo);
  p1.cprbot2(water(oo)) = p1.cprbot2(water(oo)) - delta(oo);
end

water = find(p1.ctype2 == 101 & p1.ctype == 101);
pcenter_water(water) = (p0.cprtop(water) + p0.cprbot(water))/2;
%% these are the ones we want to "raise"
oo = find(pcenter_water(water) > p1.sarta_lvl_waterOD_1(water) & p1.sarta_lvl_waterOD_1(water) > 0);
if length(oo) > 0
  delta = pcenter_water(water) - p1.sarta_lvl_waterOD_1(water);
  p1.cprtop(water(oo)) = p1.cprtop(water(oo)) - delta(oo);
  p1.cprbot(water(oo)) = p1.cprbot(water(oo)) - delta(oo);
end


