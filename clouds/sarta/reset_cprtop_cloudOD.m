function p1 = reset_cprtop_cloudOD(p0,cumsumOD,airslevels,airsheights);

%% computes cloud ODs based on formulas given by Xianglei and Xiuhong
%% see PCRTM_compute_for_AIRS_spectra.m
%% if cumsumOD < 99 then look for where cumulative cloudOD == cumsum
%%    cumsumOD > 99 then look for peak of cloud wgt fcn

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

  %% bugfix on 6.3.2013 ... this was mistakenly ptemp before early June, 2013
  Z = p2hFAST(press,airslevels,airsheights)/1000;
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

%% do the weighting functions
wgtI = zeros(size(p0.ptemp));
ii = 1;
wgtI(ii,:) = 0;
for ii = 2 : p0.nlevs(1)
  wgtI(ii,:) = wgtI(ii-1,:) + p1.sarta_lvlODice(ii,:);
end
wgtI = exp(-wgtI) .* (1-exp(-p1.sarta_lvlODice));
for ii = 1 : length(p1.stemp)
  dodo = wgtI(:,ii);
  if sum(dodo) > 0
    bop = find(dodo == max(dodo),1);
    wgtpeakIindex(ii) = bop;
    wgtpeakI_tempr(ii) = p1.ptemp(bop,ii);
    wgtpeakI(ii) = p1.plevs(bop,ii);
  else
    wgtpeakIindex(ii) = -1;
    wgtpeakI_tempr(ii) = -9999;
    wgtpeakI(ii) = -9999;
  end
end

wgtW = zeros(size(p0.ptemp));
ii = 1;
wgtI(ii,:) = 0;
for ii = 2 : p0.nlevs(1)
  wgtW(ii,:) = wgtW(ii-1,:) + p1.sarta_lvlODwater(ii,:);
end
wgtW = exp(-wgtW) .* (1-exp(-p1.sarta_lvlODwater));
for ii = 1 : length(p1.stemp)
  dodo = wgtW(:,ii);
  if sum(dodo) > 0
    bop = find(dodo == max(dodo),1);
    wgtpeakWindex(ii) = bop;
    wgtpeakW_tempr(ii) = p1.ptemp(bop,ii);
    wgtpeakW(ii) = p1.plevs(bop,ii);
  else
    wgtpeakWindex(ii) = -1;
    wgtpeakW_tempr(ii) = -9999;
    wgtpeakW(ii) = -9999;
  end
end

p1.sarta_wgtI     = wgtI;
p1.sarta_wgtW     = wgtW;

p1.sarta_index_wgtpeakW = wgtpeakWindex;
p1.sarta_index_wgtpeakI = wgtpeakIindex;

p1.sarta_wgtpeakW = wgtpeakW;
p1.sarta_wgtpeakI = wgtpeakI;

%%%%%%%%%%%%%%%%%%%%%%%%%
iJunk = -1;
if iJunk > 0
  %% this shows that wgtpeakW < sarta_lvl_waterOD_1 ie higher up!!!!!
  plot(p1.sarta_wgtpeakW,p1.sarta_lvl_waterOD_1,'o',p1.sarta_wgtpeakI,p1.sarta_lvl_iceOD_1,'rx',wgtpeakI,wgtpeakI)
   axis([0 1000 0 1000])

  plot(wgtpeakI,p1.cprtop,'o',wgtpeakW,p1.cprtop2,'ro',wgtpeakI,wgtpeakI); axis([0 1000 0 1000])
  dn = -1000 : 10 : +1000; nn = hist(wgtpeakI-p1.cprtop,dn); plot(dn,nn)
  dn = -1000 : 10 : +1000; nn = hist(wgtpeakW-p1.cprtop2,dn); plot(dn,nn)

  %% do simple RT
  ii = p0.nlevs(1);
  radI(ii,:) = ttorad(1231*ones(size(p1.stemp)),p1.stemp);
  for ii = p0.nlevs(1)-1 : -1 : +1
    tau = exp(-p1.sarta_lvlODice(ii,:));
    radI(ii,:) = radI(ii+1,:) .* tau + ttorad(1231*ones(size(p1.stemp)),p1.ptemp(ii,:)) .* (1 - tau);
  end

  ii = p0.nlevs(1);
  radW(ii,:) = ttorad(1231*ones(size(p1.stemp)),p1.stemp);
  for ii = p0.nlevs(1)-1 : -1 : +1
    tau = exp(-p1.sarta_lvlODwater(ii,:));
    radW(ii,:) = radW(ii+1,:) .* tau + ttorad(1231*ones(size(p1.stemp)),p1.ptemp(ii,:)) .* (1 - tau);
  end

  junkprof = 1 - exp(-(1:double(p0.nlevs(1)))/double(p0.nlevs(1)));  %% exponential atmosphere
  ii = p0.nlevs(1);
  radJ(ii,:) = ttorad(1231*ones(size(p1.stemp)),p1.stemp);
  for ii = p0.nlevs(1)-1 : -1 : +1
    tau = exp(-junkprof(ii))*ones(size(p1.stemp));
    radJ(ii,:) = radJ(ii+1,:) .* tau + ttorad(1231*ones(size(p1.stemp)),p1.ptemp(ii,:)) .* (1 - tau);
  end

  lala = sum(p1.sarta_lvlODice);
  scatter(ttorad(1231,wgtpeakI_tempr),radI(1,:),30,lala,'filled'); colorbar
  lala = sum(p1.sarta_lvlODwater);
  scatter(ttorad(1231,wgtpeakW_tempr),radW(1,:),30,lala,'filled'); colorbar
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now do the actual replacement of cprtop/cprbot with above info
p1 = do_the_reset_cprtop_cloudOD(p0,p1,cumsumOD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

