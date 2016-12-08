function [hx,p0ALL,p_co2_n2o_co_ch4_pcrtm] = pad_upper_atm(h,p0ALL_inputLVLS);

%% input
%%   h,p0ALL_inputLVLS are input header, profile (levels)
%%   typically p0ALL_inputLVLS geophysical is ERA or ECM, upto 1 mb
%% output
%%   hx = h, but with pmin reset to 0.005 mb
%%   p0ALL has US STd profiles for ptemp,gas_1 and gas_3 tacked on   between 1 mb to 0.005 mb
%%                                 cc,ciwc,clwc          set to zero between 1 mb to 0.005 mb
%%                                 nlevs typically increased by a few points (eg from 60 to 64 for ERA)
%%                                 txover,gxover   reset
%%      quite a few fields removed, as they will be reset while running this code eg rad_allsky_std, rad_allsky etc
%%   p_co2_n2o_co_ch4_pcrtm are the PCRTM CO2,N2O,CO,CH4 profile at the same pressure levels as are now in p0ALL, at default gnd 385.84,1.843 ppmv 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i2 = -1;  %% assume no CO2 in input profile
i4 = -1;  %% assume no N2O in input profile
i5 = -1;  %% assume no CO  in input profile
i6 = -1;  %% assume no CH4 in input profile

p0ALL = p0ALL_inputLVLS;

%% does part of the jobs of get_sarta_clear2, get_sarta_clear2 by adding in UpperAtm info (except for CO2/CH4)
%% ie almost like p_add_co2_ch4_complete, from which it is copied ... but does NOT add in CO2/CH4

pcrtm_p = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/Atm_prof/lev-101_nMol-6/pbnd.dat');
pcrtm = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/par_constant.dat');
%% should be CO2 N2O CO CH4

%% see AFGL proflies in
%% DFAFGL in /asl/packages/klayersV205/Src/incLAY.f is
%%    /asl/packages/klayersV205/Data/glatm.dat -> glatm_16Aug2010.dat
%% [h,ha,p,pa] = rtpread('/asl/packages/klayersV205/Data/adafgl_16Aug2010_ip.rtp');
load_std_profile

%% need to convert ppmv to g/g since that is what ECMWF brings
%% /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
H2O = ppmv2gg(Pressure,Temperature,H2O,18);
O3  = ppmv2gg(Pressure,Temperature,O3,48);

woop = p0ALL.plevs;
woop(woop < 0) = NaN;
minp1 = min(p0ALL.plevs(:));
minp2 = nanmin(woop(:));
if minp1 ~= minp2
  error('oops different minp???')
end

%% initial set of xover
pxx.txover      = ones(size(p0ALL.stemp)) * minp1;
pxx.txover      = ones(size(p0ALL.stemp)) * minp1;
for ix = 1 : h.ngas
  pxx.gxover(ix,:) = ones(size(p0ALL.stemp)) * minp1;
end

%% initial replacement with current profile
newpoints = find(pcrtm_p < minp1);
[mmjunk,nnjunk] = size(p0ALL.ptemp);
pxx.gas_1 = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.gas_1(length(newpoints)+1:mmjunk+length(newpoints),:) = p0ALL.gas_1;
pxx.gas_3 = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.gas_3(length(newpoints)+1:mmjunk+length(newpoints),:) = p0ALL.gas_3;
pxx.ptemp = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.ptemp(length(newpoints)+1:mmjunk+length(newpoints),:) = p0ALL.ptemp;
pxx.plevs = ones(mmjunk+length(newpoints),nnjunk) * -9999; pxx.plevs(length(newpoints)+1:mmjunk+length(newpoints),:) = p0ALL.plevs;
if isfield(p0ALL,'ciwc')
  pxx.ciwc = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.ciwc(length(newpoints)+1:mmjunk+length(newpoints),:) = p0ALL.ciwc;
  pxx.clwc = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.clwc(length(newpoints)+1:mmjunk+length(newpoints),:) = p0ALL.clwc;
  pxx.cc = ones(mmjunk+length(newpoints),nnjunk) *   0.0; pxx.cc(length(newpoints)+1:mmjunk+length(newpoints),:)   = p0ALL.cc;
end
if isfield(p0ALL,'sarta_lvlODice')
  pxx.sarta_lvlODice   = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.sarta_lvlODice(length(newpoints)+1:mmjunk+length(newpoints),:)   = p0ALL.sarta_lvlODice;
  pxx.sarta_lvlODwater = ones(mmjunk+length(newpoints),nnjunk) * 0.0; pxx.sarta_lvlODwater(length(newpoints)+1:mmjunk+length(newpoints),:) = p0ALL.sarta_lvlODwater;
end

pxx.nlevs = p0ALL.nlevs + length(newpoints);
pxx.plevs(1:length(newpoints),:) = pcrtm_p(newpoints) * ones(1,nnjunk);

%% add in WV, O3, T  profile from 1mb to TOA
pxx.txover    = ones(size(p0ALL.stemp)) * min(pcrtm_p);
pxx.gxover(1,:) = ones(size(p0ALL.stemp)) * min(pcrtm_p);
pxx.gxover(2,:) = ones(size(p0ALL.stemp)) * min(pcrtm_p);
for iijunk = 1 : length(p0ALL.stemp)
  boo = find(Pressure <= p0ALL.plevs(1,iijunk));
  boo = [boo(1)-1 boo];

  woo  = interp1(log(Pressure(boo)),Temperature(boo),log(pcrtm_p(newpoints)),'linear','extrap');
  woo0 = interp1(log(Pressure(boo)),Temperature(boo),log(p0ALL.plevs(1,iijunk)),  'linear','extrap');
  offset = p0ALL.ptemp(1,iijunk) - woo0;
  pxx.ptemp(1:length(newpoints),iijunk) = woo + offset;

  woo  = interp1(log(Pressure(boo)),H2O(boo),log(pcrtm_p(newpoints)),'linear','extrap');
  woo0 = interp1(log(Pressure(boo)),H2O(boo),log(p0ALL.plevs(1,iijunk)),  'linear','extrap');
  offset = p0ALL.gas_1(1,iijunk) ./ woo0;
  pxx.gas_1(1:length(newpoints),iijunk) = woo * offset;

  woo  = interp1(log(Pressure(boo)),O3(boo),log(pcrtm_p(newpoints)),'linear','extrap');
  woo0 = interp1(log(Pressure(boo)),O3(boo),log(p0ALL.plevs(1,iijunk)),    'linear','extrap');
  offset = p0ALL.gas_3(1,iijunk) ./ woo0;
  pxx.gas_3(1:length(newpoints),iijunk) = woo * offset;
end

p0ALL.gxover = pxx.gxover;
p0ALL.txover = pxx.txover;
p0ALL.gas_1 = pxx.gas_1;
p0ALL.gas_3 = pxx.gas_3;
p0ALL.ptemp = pxx.ptemp;
p0ALL.plevs = pxx.plevs;
p0ALL.nlevs = pxx.nlevs;
if isfield(p0ALL,'ciwc')
  p0ALL.ciwc = pxx.ciwc;
  p0ALL.clwc = pxx.clwc;
  p0ALL.cc   = pxx.cc;
end
if isfield(p0ALL,'sarta_lvlODice')
  p0ALL.sarta_lvlODice   = pxx.sarta_lvlODice;
  p0ALL.sarta_lvlODwater = pxx.sarta_lvlODwater;
end

hx      = h;
hx.pmin = min(pcrtm_p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add in CO2/N2O/CH4
%% see p_add_co2_ch4_simple.m
pxjunk = struct;
pxjunk.stemp = p0ALL.stemp;
pxjunk.nlevs = p0ALL.nlevs;
pxjunk.plevs = p0ALL.plevs;
pxjunk.ptemp = p0ALL.ptemp;

if length(intersect(hx.glist,2)) == 1
  disp('oh oh : you have supplied input CO2 profile which will be overwritten here!');
  i2 = +1;
end
if length(intersect(hx.glist,4)) == 1
  disp('oh oh : you have supplied input N2O profile which will be overwritten here!');
  i4 = +1;
end
if length(intersect(hx.glist,5)) == 1
  disp('oh oh : you have supplied input CO   profile which will be overwritten here!');
  i5 = +1;
end  
if length(intersect(hx.glist,6)) == 1
  disp('oh oh : you have supplied input CH4 profile which will be overwritten here!');
  i6 = +1;
end  

for iijunk = 1 : length(pxjunk.stemp)
  nlevs = pxjunk.nlevs(iijunk);
  boo   = pxjunk.plevs(:,iijunk);

  %% co2 = 385.848 ppm at GND
  junk = interp1(log10(pcrtm_p),pcrtm(:,1),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_2(:,iijunk) = junky;

  %% n2o
  junk = interp1(log10(pcrtm_p),pcrtm(:,2),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_4(:,iijunk) = junky;	    

  %% co
  junk = interp1(log10(pcrtm_p),pcrtm(:,3),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_5(:,iijunk) = junky;	    

  %% ch4 = 1.843 ppm at GND
  junk = interp1(log10(pcrtm_p),pcrtm(:,4),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_6(:,iijunk) = junky;	    
end
p_co2_n2o_co_ch4_pcrtm = pxjunk;
%% add in CO2/N2O/CO/CH4

if i2 > 0
  thediff = 1 - p0ALL.gas_2 ./ pxjunk.gas_2;
  thediff = nanmean(nanmean(thediff));
  fprintf(1,'mean ratio diff   1 - (input CO2/PCRTM CO2) = %8.6f \n',thediff);
  if (abs(thediff-1) > eps)
    disp('could be that PCRTM default is 385.848 ppm and you had previously asked for eg 385 ppm')
  end
end
if i2 > 0
  thediff = 1 - p0ALL.gas_4 ./ pxjunk.gas_4;
  thediff = nanmean(nanmean(thediff));
  fprintf(1,'mean ratio diff   1 - (input N2O/PCRTM N2O) = %8.6f \n',thediff);
end
if i2 > 0
  thediff = 1 - p0ALL.gas_5 ./ pxjunk.gas_5;
  thediff = nanmean(nanmean(thediff));
  fprintf(1,'mean ratio diff   1 - (input CO/PCRTM CO) = %8.6f \n',thediff);
end
if i6 > 0
  thediff = 1 - p0ALL.gas_6 ./ pxjunk.gas_6;
  thediff = nanmean(nanmean(thediff));
  fprintf(1,'mean ratio diff   1 - (input CH4/PCRTM CH4) = %8.6f \n',thediff);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
did not find field nrho
did not find field rad_allsky_std
did not find field rad_clrsky
did not find field ncol
did not find field pcrtm_co2_used
did not find field overlap
did not find field pcrtm_iceOD
did not find field pcrtm_iceDME
did not find field pcrtm_iceCTOP
did not find field pcrtm_waterOD
did not find field pcrtm_waterDME
did not find field pcrtm_waterCTOP
did not find field pcrtm_lvlODice
did not find field pcrtm_lvlODwater
did not find field pcrtm_iceODX
did not find field pcrtm_waterODX
did not find field rcalc_std
%}

if isfield(p0ALL,'rad_allsky_std')
  p0ALL = rmfield(p0ALL,'rad_allsky_std');
end
if isfield(p0ALL,'rad_clrsky')
  p0ALL = rmfield(p0ALL,'rad_clrsky');
end

if isfield(p0ALL,'pcrtm_iceOD')
  p0ALL = rmfield(p0ALL,'pcrtm_iceOD');
end
if isfield(p0ALL,'pcrtm_iceDME')
  p0ALL = rmfield(p0ALL,'pcrtm_iceDME');
end
if isfield(p0ALL,'pcrtm_iceCTOP')
  p0ALL = rmfield(p0ALL,'pcrtm_iceCTOP');
end
if isfield(p0ALL,'pcrtm_lvlODice')
  p0ALL = rmfield(p0ALL,'pcrtm_lvlODice');
end
if isfield(p0ALL,'pcrtm_iceODX')
  p0ALL = rmfield(p0ALL,'pcrtm_iceODX');
end

if isfield(p0ALL,'pcrtm_waterOD')
  p0ALL = rmfield(p0ALL,'pcrtm_waterOD');
end
if isfield(p0ALL,'pcrtm_waterDME')
  p0ALL = rmfield(p0ALL,'pcrtm_waterDME');
end
if isfield(p0ALL,'pcrtm_waterCTOP')
  p0ALL = rmfield(p0ALL,'pcrtm_waterCTOP');
end
if isfield(p0ALL,'pcrtm_lvlODwater')
  p0ALL = rmfield(p0ALL,'pcrtm_lvlODwater');
end
if isfield(p0ALL,'pcrtm_waterODX')
  p0ALL = rmfield(p0ALL,'pcrtm_waterODX');
end


