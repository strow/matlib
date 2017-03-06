function p1 = compute_cloud_OD_moleculespercm2(p0,airslevels,airsheights);

disp('ASSUMES LEVELS PROFILE ... computing level ODs and molecules/cm2')

%% basically same as reset_cprtop_cloudOD.m
%% basically same as compute_cloud_OD.m

%% computes cloud ODs based on formulas given by Xianglei and Xiuhong
%%   see PCRTM_compute_for_AIRS_spectra.m
%%   if cumsumOD < 99 then look for where cumulative cloudOD == cumsum
%%      cumsumOD > 99 then look for peak of cloud wgt fcn

%% computes molecules/cm2 based on klayers by Scott Hannon
%{
A layer amount of the type required by kCARTA can be calculated using the
equation:

   A = den_ref * (Pavg/Tavg) * MRavg * delta_z

where
   den_ref = 1.2027E-12 (kmole/cm^2)*(K/mb)*(1/ppmv)*(1/m)
      Pavg, Tavg, and delta_z are as defined earlier, and
         MRavg = the average mixing ratio.

with MRavg perhaps being approximated by:

   MRavg = (D1*MR1 + D2*MR2)/(D1 + D2)

SO eg for 1.202E-12 *1e4 kmol/m2,    273 K, 1013 mb, 1e6 ppmv (so ratio 1), 1 m
gives 0.0446 kmol/m3 or 44.6 moles/m3 or 2.67e25 molecules/m3 == loschmidt number
%}

%{
testing
suppose p0 = 92 levels profile with gas_1,plevs,ptemp,cc,ciwc,clwc
p = load('airslevels.dat');
z = load('airsheights.dat');
p1 = compute_cloud_OD_moleculespercm2(p0,p,z);
plot(abs(nansum(p1.q_gas1(1:91,:))))

%}

p1 = p0;
%p1.orig_ctop  = p1.cprtop;
%p1.orig_ctop2 = p1.cprtop2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

den_ref = 1.2027E-12; %% (kmole/cm^2)*(K/mb)*(1/ppmv)*(1/m)
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
  gas_2 = p0.gas_2(:,ii);  
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
  %% convert g/g to ppmv
  clwcMR = toppmv(press,ptemp,clwc,18.0*1.000,21);
  ciwcMR = toppmv(press,ptemp,ciwc,18.0*0.912,21);  
  wvMR   = toppmv(press,ptemp,gas_1,18.0*0.912,21);
  
  z = p2h(press);    %% in meters
  dz = abs(diff(z)); %% in meters
  clwcMR_avg = (clwcMR(1:end-1) + clwcMR(2:end))/2.0;
  ciwcMR_avg = (ciwcMR(1:end-1) + ciwcMR(2:end))/2.0;
  gas1_avg   = (wvMR(1:end-1)   + wvMR(2:end))/2.0;    
  p_avg      = (press(1:end-1) + press(2:end))/2.0;
  t_avg      = (ptemp(1:end-1) + ptemp(2:end))/2.0;  

  p1.q_ciwc(1:length(p_avg),ii) = den_ref * 1e4 * p_avg./t_avg .* ciwcMR_avg .* dz;  %% kmol/m2
  p1.q_clwc(1:length(p_avg),ii) = den_ref * 1e4 * p_avg./t_avg .* clwcMR_avg .* dz;  %% kmol/m2
  p1.q_gas1(1:length(p_avg),ii) = den_ref * 1e4 * p_avg./t_avg .* gas1_avg   .* dz;  %% kmol/m2
  
  p1.q_ciwc(length(p_avg)+1,ii) = -9999;
  p1.q_clwc(length(p_avg)+1,ii) = -9999;
  p1.q_gas1(length(p_avg)+1,ii) = -9999;

  p1.q_ciwc(:,ii) =  p1.q_ciwc(:,ii) / 1e4 * 6.023e23 * 1000;  %% convert to kmol/cm2  and then molecules
  p1.q_clwc(:,ii) =  p1.q_clwc(:,ii) / 1e4 * 6.023e23 * 1000;  %% convert to kmol/cm2  and then molecules
  p1.q_gas1(:,ii) =  p1.q_gas1(:,ii) / 1e4 * 6.023e23 * 1000;  %% convert to kmol/cm2  and then molecules  
  
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mmwater_rtp
% Conversion factor (equivalent mm of liquid water per molecules/cm^2)
% cfac = 1 (cm^3/gram) * ~18/Navadagro (grams/molecule) * 10 (mm/cm)

% Avagadro's number (molecules per mole), mass of water per molecule (AMU)
navagadro=6.02214199E+23;
mass=18.015;
cfac=10*mass/navagadro;

[mm,nn] = size(p1.plevs);
plot(abs(nansum(p1.q_gas1(1:mm-1,:))));
p1.mmw = cfac * nansum(p1.q_gas1(1:mm-1,:));

%% each water/ice molecule = 18g/avog number
p1.totalmass_ciwc = nansum(p1.q_ciwc(1:mm-1,:)) * 18/navagadro * 1e4;   %% amt in molecules/cm2 --> mass in g/cm2 --> g/m2
p1.totalmass_clwc = nansum(p1.q_clwc(1:mm-1,:)) * 18/navagadro * 1e4;   %% amt im molecules/cm2 --> mass in g/cm2 --> g/m2

plot(p1.mass_ciwc,p1.cngwat,'b.',p1.mass_clwc,p1.cngwat2,'r.')     %% should be roughly the same
  xlabel('mass ciwc, mass clwc (g/m2)'); ylabel('cngwat(ice/b) cngwat2(water/r)');
line([0 500],[0 500],'color','k','linewidth',3)
axis([0 500 0 500])

%% this is conversion of OD to g/m2
plot(p1.sarta_lvlODice(:,1:100) , p1.q_ciwc(:,1:100)*18/6e23*1e4,'.'); axis([0 0.25 0 0.25])
plot(p1.sarta_lvlODice(:,1:100) ./ p1.q_ciwc(:,1:100)*18/6e23*1e4,'.'); axis([0 0.25 0 0.25])

zooA = p1.sarta_lvlODice(1:mm-1,:);
zooB = p1.q_ciwc(1:mm-1,:)*18/6e23*1e4;
zooA = zooA(:);
zooB = zooB(:);
ix = find(zooA > 0 & zooA > 0);
semilogy(0:0.01:100,hist(zooA(ix) ./ zooB(ix),0:0.01:100))