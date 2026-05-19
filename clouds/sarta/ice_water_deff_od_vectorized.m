function [p1,iceOD,waterOD,extrajunk] = ice_water_deff_od(p0,airslevels,airslayers,airsheights,iNew_or_Orig_CXWC2OD)

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

ciwc  = p0.ciwc;
clwc  = p0.clwc;
cc    = p0.cc;
ptemp = p0.ptemp;
gas_1 = p0.gas_1;
press = p0.plevs;

pressN = press(:,1);
%% code works assming press(1) < press(2) etc ie TOA is at LVL1, GND = LVLN
if pressN(1) > pressN(2)
  clwc  = flipud(clwc);
  ciwc  = flipud(ciwc);
  cc    = flipud(cc);
  ptemp = flipud(ptemp);
  gas_1 = flipud(gas_1);
  press = flipud(press);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
tcld = ptemp - 273.16;                    %% change to deg C
boo  = find(tcld < -50); tcld(boo) = -50; %% set minimum as -50 C
boo  = find(tcld > -25); tcld(boo) = -25; %% set maximum as -50 C

[mm,nn] = size(ptemp);
tav    = (ptemp(1:mm-1,:)+ptemp(2:mm,:))/2; tav(mm,:) = ptemp(mm,:);
pav    = (press(1:mm-1,:)+press(2:mm,:))/2; pav(mm,:) = press(mm,:);
scaleH = R*tav/g/1000;    %% in km
dz     = scaleH(1:mm-1,:) .* log(press(2:mm,:)./press(1:mm-1,:)); dz(mm,:) = 0;
Z      = cumsum(dz,1);
diffZ  = abs(diff(Z,1)); diffZ(mm,:) = diffZ(mm-1,:);

%% bugfix on 6.3.2013 ... this was mistakenly ptemp before early June, 2013
Z     = p2hFAST(press,airslevels,airslayers,airsheights)/1000;
%whos press airslevels airsheights
diffZ = abs(diff(Z,1)); diffZ(mm,:) = diffZ(mm-1,:);

if nargout == 4
  extrajunk.tav = tav;
  extrajunk.pav = pav;
  extrajunk.scaleH = scaleH;
  extrajunk.Z = Z;
  extrajunk.diffZ = diffZ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
p1.sarta_lvlZ = Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cld_de_ice = c0 + c1 * tcld + c2 * tcld.^2 + c3 * tcld.^3;

% compute ice cloud optical depth from Ebert and Curry (1992, J. Geophys. Res.) 
qi = ciwc ./ ptemp .* press *100/R * 1e3;  %%change IWC from kg/kg to g/m^3
if iNew_or_Orig_CXWC2OD == -1
  iceOD = (0.003448 + 2.431./cld_de_ice) .* qi ./ cc .* diffZ *1e3;   %% ORIG
elseif iNew_or_Orig_CXWC2OD == 0
  iceOD = (0.003448 + 2.431./cld_de_ice) .* qi ./ cc .* diffZ *1e3;   %% XIUHONG
  iceOD(cc < 1e-3) = 0.0;
elseif iNew_or_Orig_CXWC2OD == +1
  iceOD = (0.003448 + 2.431./cld_de_ice) .* qi .* cc .* diffZ *1e3;   %% SERGIO
end

bad = find(ciwc < theeps);
iceOD(bad) = 0.0;

%figure(3); plot(cld_de_ice); title('cld_de ice'); disp('ret'); pause
%figure(3); plot(qi);        title('qi'); disp('ret'); pause
%figure(3); plot(diffZ);     title('diffZ'); disp('ret'); pause
%figure(3); plot(cc);        title('cc'); disp('ret'); pause
%figure(3); plot(iceOD);     title('iceOD 0'); disp('ret'); pause

bad = find(isnan(cld_de_ice)); cld_de_ice(bad) = -9999;
p1.sarta_lvlDMEice = cld_de_ice;

bad = find(isnan(iceOD) | isinf(iceOD)); iceOD(bad) = 0;
p1.sarta_lvlODice = iceOD;

%figure(3); plot(iceOD);     title('iceOD F'); disp('ret'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cld_de_liq = 20 * ones(size(press));

% compute water cloud optical depth, % ECMWF technical report, Klein S. A., 1999
qw = clwc ./ ptemp .* press *100/R *1e3;  %change liquid water content from kg/kg to g/m^3
if iNew_or_Orig_CXWC2OD == -1
  waterOD = 3 * qw ./ cld_de_liq ./cc .*diffZ  *1e3;  %% ORIG
elseif iNew_or_Orig_CXWC2OD == 0  
  waterOD = 3 * qw ./ cld_de_liq ./cc .*diffZ  *1e3;  %% XIUHONG
  waterOD(cc < 1e-3) = 0.0;
elseif iNew_or_Orig_CXWC2OD == +1  
  waterOD = 3 * qw ./ cld_de_liq .*cc .*diffZ  *1e3;  %% SERGIO
end

bad = find(isnan(cld_de_liq)); cld_de_liq(bad) = -9999;
p1.sarta_lvlDMEwater = cld_de_liq;

bad = find(isnan(waterOD) | isinf(waterOD)); waterOD(bad) = 0;
p1.sarta_lvlODwater = waterOD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
