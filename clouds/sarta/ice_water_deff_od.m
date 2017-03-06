function [p1,iceOD,waterOD] = ice_water_deff_od(p0,airslevels,airsheights,ii,iNew_or_Orig_CXWC2OD)

if nargin == 4
  iNew_or_Orig_CXWC2OD =  0;  %%% change to OD = blah * qBlah / cc * diffZ; OD(cc < 1e-3) = 0 WHAT PCRTM DOES
  iNew_or_Orig_CXWC2OD = +1;  %%% change to OD = blah * qBlah * cc * diffZ                    Mar 2017 SERGIO
  iNew_or_Orig_CXWC2OD = -1;  %%% stick  to OD = blah * qBlah / cc * diffZ                    Pre March 2017  DEFAULT
end

theeps = 1e-15;

R = 287.05;
g = 9.806;

% coefficents from S-C Ou, K-N. Liou, Atmospheric Research
% 35(1995):127-138.
% for computing ice cloud effective size
c0 = 326.3;
c1 = 12.42;
c2 = 0.197;
c3 = 0.0012;

p1 = p0;

ciwc = p0.ciwc(:,ii);
clwc = p0.clwc(:,ii);
cc   = p0.cc(:,ii);
ptemp = p0.ptemp(:,ii);
gas_1 = p0.gas_1(:,ii);
press = p0.plevs(:,ii);

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

tav = (ptemp(1:end-1)+ptemp(2:end))/2; tav(end) = ptemp(end);
pav = (press(1:end-1)+press(2:end))/2; pav(end) = press(end);
scaleH = R*tav/g/1000;    %% in km
dz = scaleH .* log(press(2:end)./press(1:end-1)); dz(length(dz)+1) = 0;
Z = cumsum(dz);
diffZ = abs(diff(Z)); diffZ(length(diffZ)+1) = diffZ(length(diffZ));

%% bugfix on 6.3.2013 ... this was mistakenly ptemp before early June, 2013
Z = p2hFAST(press,airslevels,airsheights)/1000;
diffZ = abs(diff(Z)); diffZ(length(diffZ)+1) = diffZ(length(diffZ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p1.sarta_lvlZ(:,ii) = Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cldde_ice = c0 + c1 * tcld + c2 * tcld.^2 + c3 * tcld.^3;

% compute ice cloud optical depth from Ebert and Curry (1992, J. Geophys. Res.) 
qi = ciwc ./ ptemp .* press *100/R * 1e3;  %%change IWC from kg/kg to g/m^3
if iNew_or_Orig_CXWC2OD == -1
  iceOD = (0.003448 + 2.431./cldde_ice) .* qi ./ cc .* diffZ *1e3;   %% ORIG
elseif iNew_or_Orig_CXWC2OD == 0
  iceOD = (0.003448 + 2.431./cldde_ice) .* qi ./ cc .* diffZ *1e3;   %% XIUHONG
  iceOD(cc < 1e-3) = 0.0;
elseif iNew_or_Orig_CXWC2OD == +1
  iceOD = (0.003448 + 2.431./cldde_ice) .* qi .* cc .* diffZ *1e3;   %% SERGIO
end

bad = find(ciwc < theeps);
iceOD(bad) = 0.0;

%figure(3); plot(cldde_ice); title('cldde ice'); disp('ret'); pause
%figure(3); plot(qi);        title('qi'); disp('ret'); pause
%figure(3); plot(diffZ);     title('diffZ'); disp('ret'); pause
%figure(3); plot(cc);        title('cc'); disp('ret'); pause
%figure(3); plot(iceOD);     title('iceOD 0'); disp('ret'); pause

bad = find(isnan(cldde_ice)); cldde_ice(bad) = -9999;
p1.sarta_lvlDMEice(:,ii) = cldde_ice;

bad = find(isnan(iceOD) | isinf(iceOD)); iceOD(bad) = 0;
p1.sarta_lvlODice(:,ii) = iceOD;

%figure(3); plot(iceOD);     title('iceOD F'); disp('ret'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cldde_liq = 20 * ones(size(press));

% compute water cloud optical depth, % ECMWF technical report, Klein S. A., 1999
qw = clwc ./ ptemp .* press *100/R *1e3;  %change liquid water content from kg/kg to g/m^3
if iNew_or_Orig_CXWC2OD == -1
  waterOD = 3 * qw ./ cldde_liq ./cc .*diffZ  *1e3;  %% ORIG
elseif iNew_or_Orig_CXWC2OD == 0  
  waterOD = 3 * qw ./ cldde_liq ./cc .*diffZ  *1e3;  %% XIUHONG
  waterOD(cc < 1e-3) = 0.0;
elseif iNew_or_Orig_CXWC2OD == +1  
  waterOD = 3 * qw ./ cldde_liq .*cc .*diffZ  *1e3;  %% SERGIO
end

bad = find(isnan(cldde_liq)); cldde_liq(bad) = -9999;
p1.sarta_lvlDMEwater(:,ii) = cldde_liq;

bad = find(isnan(waterOD) | isinf(waterOD)); waterOD(bad) = 0;
p1.sarta_lvlODwater(:,ii) = waterOD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
