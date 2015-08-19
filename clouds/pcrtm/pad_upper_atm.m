function [hx,p0ALL] = pad_upper_atm(h,p0ALL_inputLVLS);

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
pxx.gxover(1,:) = ones(size(p0ALL.stemp)) * minp1;
pxx.gxover(2,:) = ones(size(p0ALL.stemp)) * minp1;

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


