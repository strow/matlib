%% take care of ALL land/ocean
if iDayNight == 0 & iLandOcean == 0
  ix = find(p.robs1(oo,:) > 0);
elseif iDayNight == -1 & iLandOcean == 0
  ix = find(p.robs1(oo,:) > 0 & p.solzen > 90);
elseif iDayNight == +1 & iLandOcean == 0
  ix = find(p.robs1(oo,:) > 0 & p.solzen <= 90);

%% take care of ALL LAND
elseif iDayNight == 0 & iLandOcean == -1
  ix = find(p.robs1(oo,:) > 0 & p.landfrac > 0.99);
elseif iDayNight == -1 & iLandOcean == -1
  ix = find(p.robs1(oo,:) > 0 & p.solzen > 90 & p.landfrac > 0.99);
elseif iDayNight == +1 & iLandOcean == -1
  ix = find(p.robs1(oo,:) > 0 & p.solzen <= 90 & p.landfrac > 0.99);

%% take care of ALL OCEAN
elseif iDayNight == 0 & iLandOcean == -2
  ix = find(p.robs1(oo,:) > 0 & p.landfrac < 0.01);
elseif iDayNight == -1 & iLandOcean == -2
  ix = find(p.robs1(oo,:) > 0 & p.solzen > 90 & p.landfrac < 0.01);
elseif iDayNight == +1 & iLandOcean == -2
  ix = find(p.robs1(oo,:) > 0 & p.solzen <= 90 & p.landfrac < 0.01);

else
  site = transcom_match(p.rlat, p.rlon);
  if iDayNight == 0 
    ix = find(p.robs1(oo,:) > 0 & site == iLandOcean);
  elseif iDayNight == -1
    ix = find(p.robs1(oo,:) > 0 & site == iLandOcean & p.solzen > 90);
  elseif iDayNight == +1
    ix = find(p.robs1(oo,:) > 0 & site == iLandOcean & p.solzen <= 90);
  end
end

