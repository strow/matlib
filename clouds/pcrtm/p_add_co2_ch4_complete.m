if iDoKlayersHard < 0   %% simple, just add in CO2 profiles
  disp(' adding in ONLY co2/ch4 profile using current plevs')
  p_add_co2_ch4_simple   %% loads in pcrtm CO2/CH4 and pressure profiles, and interpolates them onto p.plevs

elseif iDoKlayersHard > 0   %% just add in CO2 profiles from Xiahoong plus WV,O3,T form USSTD
  disp(' adding in co2/ch4 profile using current plevs PLUS adding on info from top-of-input-plevs to 0.005 mb')
  
  p_add_co2_ch4_simple   %% loads in pcrtm CO2/CH4 and pressure profiles, and interpolates them onto p.plevs

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%% NEW %%%
  %% DFAFGL in /asl/packages/klayersV205/Src/incLAY.f is
  %%    /asl/packages/klayersV205/Data/glatm.dat -> glatm_16Aug2010.dat
  %% [h,ha,p,pa] = rtpread('/asl/packages/klayersV205/Data/adafgl_16Aug2010_ip.rtp');
  load_std_profile
  	  
  %% need to convert ppmv to g/g since that is what ECMWF brings
  %% /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
  H2O = ppmv2gg(Pressure,Temperature,H2O,18);
  O3  = ppmv2gg(Pressure,Temperature,O3,48);

  woop = pxjunk.plevs;
  woop(woop < 0) = NaN;
  minp1 = min(pxjunk.plevs(:));
  minp2 = nanmin(woop(:));
  if minp1 ~= minp2
    error('oops different minp???')
  end

  clear pxx
  
  %% initial set of xover
  pxx.txover      = ones(size(pxjunk.stemp)) * minp1;
  pxx.txover      = ones(size(pxjunk.stemp)) * minp1;
  pxx.gxover(1,:) = ones(size(pxjunk.stemp)) * minp1;
  pxx.gxover(2,:) = ones(size(pxjunk.stemp)) * min(pcrtm_p);  
  pxx.gxover(3,:) = ones(size(pxjunk.stemp)) * minp1;
  pxx.gxover(4,:) = ones(size(pxjunk.stemp)) * min(pcrtm_p);  

  %% initial replacement with current profile
  newpoints = find(pcrtm_p < minp1);
  [mmjunk,nnjunk] = size(pxjunk.ptemp);
  pxx.gas_1 = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.gas_1(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.gas_1;
  pxx.gas_2 = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.gas_2(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.gas_2;
  pxx.gas_3 = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.gas_3(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.gas_3;
  pxx.gas_6 = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.gas_6(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.gas_6;  
  pxx.ptemp = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.ptemp(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.ptemp;
  pxx.plevs = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.plevs(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.plevs;  
  if isfield(pxjunk,'ciwc')
    pxx.ciwc = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.ciwc(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.ciwc;
    pxx.clwc = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.clwc(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.clwc;
    pxx.cc = ones(mmjunk+length(newpoints),nnjunk) *   0.0; pxx.cc(length(newpoints)+1:mmjunk+length(newpoints),:)   = pxjunk.cc;
  end
  if isfield(p0ALL,'sarta_lvlODice')
    pxx.sarta_lvlODice   = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.sarta_lvlODice(length(newpoints)+1:mmjunk+length(newpoints),:)   = pxjunk.sarta_lvlODice;
    pxx.sarta_lvlODwater = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.sarta_lvlODwater(length(newpoints)+1:mmjunk+length(newpoints),:) = pxjunk.sarta_lvlODwater;
  end

  %% add in CO2/CH4 profile from 1mb to TOA
  pxx.nlevs = pxjunk.nlevs + length(newpoints);
  pxx.plevs(1:length(newpoints),:) = pcrtm_p(newpoints) * ones(1,nnjunk);
  pxx.gas_2(1:length(newpoints),:) = pcrtm(newpoints,1) * ones(1,nnjunk) * sarta_gas_2_6.co2/385.848;
  pxx.gas_6(1:length(newpoints),:) = pcrtm(newpoints,4) * ones(1,nnjunk) * sarta_gas_2_6.ch4/1.843;

  %% add in WV, O3, T  profile from 1mb to TOA
  pxx.txover      = ones(size(pxjunk.stemp)) * min(pcrtm_p);
  pxx.gxover(1,:) = ones(size(pxjunk.stemp)) * min(pcrtm_p);
  pxx.gxover(3,:) = ones(size(pxjunk.stemp)) * min(pcrtm_p);
  for iijunk = 1 : length(pxjunk.stemp)
    boo = find(Pressure <= pxjunk.plevs(1,iijunk));
    boo = [boo(1)-1 boo];

    woo  = interp1(log(Pressure(boo)),Temperature(boo),log(pcrtm_p(newpoints)),'linear','extrap');
    woo0 = interp1(log(Pressure(boo)),Temperature(boo),log(pxjunk.plevs(1,iijunk)),    'linear','extrap');
    offset = pxjunk.ptemp(1,iijunk) - woo0;
    pxx.ptemp(1:length(newpoints),iijunk) = woo + offset;

    woo  = interp1(log(Pressure(boo)),H2O(boo),log(pcrtm_p(newpoints)),'linear','extrap');
    woo0 = interp1(log(Pressure(boo)),H2O(boo),log(pxjunk.plevs(1,iijunk)),    'linear','extrap');
    offset = pxjunk.gas_1(1,iijunk) ./ woo0;
    pxx.gas_1(1:length(newpoints),iijunk) = woo * offset;

    woo  = interp1(log(Pressure(boo)),O3(boo),log(pcrtm_p(newpoints)),'linear','extrap');
    woo0 = interp1(log(Pressure(boo)),O3(boo),log(pxjunk.plevs(1,iijunk)),    'linear','extrap');
    offset = pxjunk.gas_3(1,iijunk) ./ woo0;
    pxx.gas_3(1:length(newpoints),iijunk) = woo * offset;

  end

  pxjunk.gxover = pxx.gxover;
  pxjunk.txover = pxx.txover;
  pxjunk.gas_1 = pxx.gas_1;
  pxjunk.gas_2 = pxx.gas_2;
  pxjunk.gas_3 = pxx.gas_3;
  pxjunk.gas_6 = pxx.gas_6;  
  pxjunk.ptemp = pxx.ptemp;
  pxjunk.plevs = pxx.plevs;
  pxjunk.nlevs = pxx.nlevs;
  if isfield(pxjunk,'ciwc')
    pxjunk.ciwc = pxx.ciwc;
    pxjunk.clwc = pxx.clwc;
    pxjunk.cc   = pxx.cc;
  end
  if isfield(p0ALL,'sarta_lvlODice')
    pxjunk.sarta_lvlODice   = pxx.sarta_lvlODice;
    pxjunk.sarta_lvlODwater = pxx.sarta_lvlODwater;
  end
  hx.pmin = min(pcrtm_p);
  %%%%%%%%%%%%%%%%%%%%%%%%%

end
