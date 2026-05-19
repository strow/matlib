function p1 = cloudOD_for100layer(p0,cumsumOD,airslevels,airslayers,airsheights);

%% computes cloud ODs based on formulas given by Xianglei and Xiuhong
%% see PCRTM_compute_for_AIRS_spectra.m
%%
%% copied from reset_cprtop_cloudOD.m, but skips the weighting fcns ....

p1 = p0;

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

  %% bugfix on 6.3.2013 ... this was mistakenly ptemp before early June, 2013
  Z = p2hFAST(press,airslevels,airslayers,airsheights)/1000;
  diffZ = abs(diff(Z)); diffZ(length(diffZ)+1) = diffZ(length(diffZ));

  % compute ice cloud optical depth from Ebert and Curry (1992, J. Geophys. Res.,  
  qi = ciwc ./ ptemp .* press *100/R * 1e3;  %%change IWC from kg/kg to g/m3  
  iceOD = (0.003448 + 2.431./cldde_ice) .*qi ./ cc .* diffZ *1e3;

  bad = find(isnan(cldde_ice)); cldde_ice(bad) = -9999;
  p1.sarta_lvlDMEice(:,ii) = cldde_ice;
  bad = find(isnan(iceOD) | isinf(iceOD)); iceOD(bad) = 0;

  %% now check things are ok ie where there is iwc, there is finite cc
  
  p1.sarta_lvlODice(:,ii) = iceOD;
  p1.sarta_lvlZ(:,ii) = Z;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  cldde_liq = 20 * ones(size(press));

  qw = clwc ./ ptemp .* press *100/R *1e3;  %change liquid water content from kg/kg to g/m^3
  % ECMWF technical report, Klein S. A., 1999
  waterOD = 3 * qw ./ cldde_liq ./cc .*diffZ  *1e3;

  bad = find(isnan(cldde_liq)); cldde_liq(bad) = -9999;
  p1.sarta_lvlDMEwater(:,ii) = cldde_liq;
  bad = find(isnan(waterOD) | isinf(waterOD)); waterOD(bad) = 0;
  p1.sarta_lvlODwater(:,ii) = waterOD;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  icesum = cumsum(iceOD);
  watersum = cumsum(waterOD);

  od_ice1_top = find(icesum >= cumsumOD,1);
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

  od_water1_top = find(watersum >= cumsumOD,1);
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

  p1.sarta_iceOD_warn(ii) = -1;  %% assume things ok between ciwc and cc
  if sum(iceOD) < eps & sum(ciwc)*1e5 > 0
    p1.sarta_iceOD_warn(ii) = +1;
  end
  p1.sarta_waterOD_warn(ii) = -1;  %% assume things ok between clwc and cc
  if sum(waterOD) < eps & sum(clwc)*1e5 > 0
    p1.sarta_waterOD_warn(ii) = +1;
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
