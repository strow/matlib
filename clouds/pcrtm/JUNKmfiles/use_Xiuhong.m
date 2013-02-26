clear P WCT ICT cc TT q o3 Ps Ts sfctype efreq emos zen_anf co2

clear p h 
disp('cleared p and h and using profile given to Xianglei in early Nov')
load xianglei_2012_05_01_00hrs

  P   = double(p.plevs);  %% 1 mb = 1 hPa

  WCT = double(p.clwc);   %% cloud liquid water content in kg/kg
  ICT = double(p.ciwc);   %% cloud ice    water content in kg/kg
  cc  = double(p.cc);     %% cloud fraction

  TT = double(p.ptemp);   %% tempertaure profile
  q  = double(p.gas_1);   %% wv profile in g/kg
  o3 = double(p.gas_3);   %% o3 profile in g/kg

  Ps = double(p.spres);   %% surface pressure in mb
  Ts = double(p.stemp);   %% surface temp in K
  
  sfctype = ones(size(p.stemp)) * -9999;   %% so that we can specify emissivity
  efreq   = double(p.efreq);
  emis    = double(p.emis);

  zen_ang = double(p.scanang);
  co2     = ones(size(p.stemp)) * (370 + (yy-2002)*2.2);

