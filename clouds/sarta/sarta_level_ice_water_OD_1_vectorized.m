function p1 = sarta_level_ice_water_OD_1_vectorized(iceOD,waterOD,cumsumOD,p0);

%% find where iceOD = cumsumOD, waterOD = cumsumOD
%% if abs(cumsumOD) < 99   then look for where cumulative cloudOD == cumsum
%%    abs(cumsumOD) = 9999 then look for peak of cloud wgt fcn, if abs(cumsumOD) == +9999 do    
%%                         STROW PICK (ice cloud can be between 0 < p < 1000 mb),
%%                                                              if abs(cumsumOD) == -9999 do 
%%                         MODIFIED STROW PICK (ice cloud can be between 0 < p < 400 mb),

p1 = p0;

ciwc  = p1.ciwc;
clwc  = p1.clwc;
cc    = p1.cc;
ptemp = p1.ptemp;
gas_1 = p1.gas_1;
press = p1.plevs;

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

icesum   = cumsum(iceOD,1);
watersum = cumsum(waterOD,1);

if abs(cumsumOD) <= 99
  xcumsumOD = abs(cumsumOD);   %% looking for the cumulative OD to be below a certain value
elseif abs(cumsumOD) == 9999
  xcumsumOD = 1.0;             %% looking for this OD value as that is typically where peak of weighting function is
end

[mm,nn] = size(ptemp);

%addpath /home/sergio/MATLABCODE
%keyboard_nowindow

iNewOrOld = +1;
if iNewOrOld > 0 
  raaInd       = (1:mm)' * ones(1,nn);
  raaXcumsumOD = xcumsumOD * ones(size(press));
  od_ice1_top  = zeros(size(press));
  od_ice1_top  = (icesum >= xcumsumOD).*raaInd;  od_ice1_top(od_ice1_top == 0) = (mm+1);
  od_ice1_top  = min(od_ice1_top,[],1);
  yes = find(od_ice1_top <= mm);
  no  = find(od_ice1_top == mm+1);
  for ii = 1 : length(yes)
    sarta_lvl_iceOD_1x(yes(ii)) = press(od_ice1_top(yes(ii)),yes(ii));
  end
  if length(no) > 0
    od_ice1_top  = zeros(size(press));
    od_ice1_top  = (iceOD > 0).*raaInd;  od_ice1_top(od_ice1_top == 0) = -9999;
    od_ice1_top  = max(od_ice1_top,[],1);
    sarta_lvl_iceOD_1x(no) = -9999;
    noo = find(od_ice1_top(no) > 0);
    for ii = 1 : length(noo)
      sarta_lvl_iceOD_1x(no(noo(ii))) = press(od_ice1_top(no(noo(ii))),no(noo(ii)));
    end
  end

  raaInd       = (1:mm)' * ones(1,nn);
  raaXcumsumOD = xcumsumOD * ones(size(press));
  od_water1_top  = zeros(size(press));
  od_water1_top  = (watersum >= xcumsumOD).*raaInd;  od_water1_top(od_water1_top == 0) = (mm+1);
  od_water1_top  = min(od_water1_top,[],1);
  yes = find(od_water1_top <= mm);
  no  = find(od_water1_top == mm+1);
  for ii = 1 : length(yes)
    sarta_lvl_waterOD_1x(yes(ii)) = press(od_water1_top(yes(ii)),yes(ii));
  end
  if length(no) > 0
    od_water1_top  = zeros(size(press));
    od_water1_top  = (waterOD > 0).*raaInd;  od_water1_top(od_water1_top == 0) = -9999;
    od_water1_top  = max(od_water1_top,[],1);
    sarta_lvl_waterOD_1x(no) = -9999;
    noo = find(od_water1_top(no) > 0);
    for ii = 1 : length(noo)
      sarta_lvl_waterOD_1x(no(noo(ii))) = press(od_water1_top(no(noo(ii))),no(noo(ii)));
    end
  end

  p1.sarta_lvl_iceOD_1   = sarta_lvl_iceOD_1x;
  p1.sarta_lvl_waterOD_1 = sarta_lvl_waterOD_1x;

else
  for ii = 1 : length(p1.stemp)
    od_ice1_top = find(icesum(:,ii) >= xcumsumOD,1);
    if length(od_ice1_top) > 0
      p1.sarta_lvl_iceOD_1(ii) = press(od_ice1_top,ii);
    else
      oo = find(iceOD(:,ii) > 0);
      if length(oo) > 0
        %% cumulative OD < 1, put cloud as low as possible
        oo = oo(end);
        od_ice1_top = oo;
        p1.sarta_lvl_iceOD_1(ii) = press(oo,ii);
      else
        od_ice1_top = NaN;
        p1.sarta_lvl_iceOD_1(ii) = -9999;
      end
    end
    
    od_water1_top = find(watersum(:,ii) >= xcumsumOD,1);
    if length(od_water1_top) > 0
      p1.sarta_lvl_waterOD_1(ii) = press(od_water1_top,ii);
    else
      oo = find(waterOD(:,ii) > 0);
      if length(oo) > 0
        %% cumulative OD < 1, put as low as possible
        oo = oo(end);
        od_water1_top = oo;
        p1.sarta_lvl_waterOD_1(ii) = press(oo,ii);
      else
        od_water1_top = NaN;
        p1.sarta_lvl_waterOD_1(ii) = -9999;
      end
    end
  end
end

p1.sarta_iceOD_warn = -1 * ones(size(p1.stemp));  %% assume things ok between ciwc and cc
sum1 = sum(iceOD,1); sum2 = sum(ciwc,1);
bad = find(sum1 < eps & sum2*1e5 > 0);
p1.sarta_iceOD_warn(bad) = +1;

p1.sarta_waterOD_warn = -1 * ones(size(p1.stemp));  %% assume things ok between clwc and cc
sum1 = sum(waterOD,1); sum2 = sum(clwc,1);
bad = find(sum1 < eps & sum2*1e5 > 0);
p1.sarta_waterOD_warn(bad) = +1;
