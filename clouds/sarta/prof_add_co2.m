prof.co2ppm = -9999 * ones(size(prof.stemp));
if run_sarta.co2ppm == -1
  % 12784 * 86400 + 27 = 1.1045e+09;
  if nanmean(prof.rtime) > 1e9
    %% /asl/matlab2012/airs/readers/xreadl1b_all.m
    [yy,mm,dd,hh] = tai2utc(profL.rtime - (12784 * 86400 + 27));
  else
    [yy,mm,dd,hh] = tai2utc(prof.rtime);
  end
  %co2   = ones(size(prof.stemp)) .* (370 + (yy-2002)*2.2);
  deltaT = (yy-2002) + (mm-1)/12 + dd/30/12;
  co2    = ones(size(prof.stemp)) .* (370 + deltaT*2.2);
elseif run_sarta.co2ppm == 0
  co2     = ones(size(prof.stemp)) * 385;
elseif run_sarta.co2ppm > 0
  co2     = ones(size(prof.stemp)) * run_sarta.co2ppm;
end

if ~isfield(p,'gas_2') & length(intersect(h.glist,2)) == 0
  fprintf(1,'adding in prof.co2ppm from run_sarta.co2ppm  = %2i \n',run_sarta.co2ppm)
  prof.co2ppm = co2;
else
  disp('in driver_sarta_cloud_rtp.m we ALREADY have gas2 in h,p')
  disp('will use this profile instead of p.co2ppm')
  prof.co2ppm = -9999 * ones(size(prof.stemp));
end
