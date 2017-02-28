function p1 = sarta_level_ice_water_OD_1(iceOD,waterOD,cumsumOD,p1x,ii);

%% find where iceOD = cumsumOD, waterOD = cumsumOD
%% if abs(cumsumOD) < 99   then look for where cumulative cloudOD == cumsum
%%    abs(cumsumOD) = 9999 then look for peak of cloud wgt fcn, if abs(cumsumOD) == +9999 do          STROW PICK (ice cloud can be between 0 < p < 1000 mb),
%%                                                              if abs(cumsumOD) == -9999 do MODIFIED STROW PICK (ice cloud can be between 0 < p < 400 mb),

p1 = p1x;

ciwc = p1.ciwc(:,ii);
clwc = p1.clwc(:,ii);
cc   = p1.cc(:,ii);
ptemp = p1.ptemp(:,ii);
gas_1 = p1.gas_1(:,ii);
press = p1.plevs(:,ii);

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

icesum = cumsum(iceOD);
watersum = cumsum(waterOD);

if abs(cumsumOD) <= 99
  xcumsumOD = abs(cumsumOD);   %% looking for the cumulative OD to be below a certain value
elseif abs(cumsumOD) == 9999
  xcumsumOD = 1.0;             %% looking for this OD value as that is typically where peak of weighting function is
end

od_ice1_top = find(icesum >= xcumsumOD,1);
if length(od_ice1_top) > 0
  p1.sarta_lvl_iceOD_1(ii) = press(od_ice1_top);
else
  oo = find(iceOD > 0);
  if length(oo) > 0
    %% cumulative OD < 1, put cloud as low as possible
    oo = oo(end);
    od_ice1_top = oo;
    p1.sarta_lvl_iceOD_1(ii) = press(oo);
  else
    od_ice1_top = NaN;
    p1.sarta_lvl_iceOD_1(ii) = -9999;
  end
end

od_water1_top = find(watersum >= xcumsumOD,1);
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

