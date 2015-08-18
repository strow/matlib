function [hx,p0ALL] = pad_upper_atm(h,p0ALL_inputLVLS);

p0ALL = p0ALL_inputLVLS;

%% does part of the jobs of get_sarta_clear2, get_sarta_clear2 by adding in UpperAtm info (except for CO2/CH4)
%% ie almost like p_add_co2_ch4_complete, from which it is copied ... but does NOT add in CO2/CH4

pcrtm_p = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/Atm_prof/lev-101_nMol-6/pbnd.dat');
pcrtm = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/par_constant.dat');
%% should be CO2 N2O CO CH4

Pressure = [...
  1.013E+03, 9.040E+02, 8.050E+02, 7.150E+02, 6.330E+02,...
  5.590E+02, 4.920E+02, 4.320E+02, 3.780E+02, 3.290E+02,...
  2.860E+02, 2.470E+02, 2.130E+02, 1.820E+02, 1.560E+02,...
  1.320E+02, 1.110E+02, 9.370E+01, 7.890E+01, 6.660E+01,...
  5.650E+01, 4.800E+01, 4.090E+01, 3.500E+01, 3.000E+01,...
  2.570E+01, 1.763E+01, 1.220E+01, 8.520E+00, 6.000E+00,...
  4.260E+00, 3.050E+00, 2.200E+00, 1.590E+00, 1.160E+00,...
  8.540E-01, 4.560E-01, 2.390E-01, 1.210E-01, 5.800E-02,...
  2.600E-02, 1.100E-02, 4.400E-03, 1.720E-03, 6.880E-04,...
  2.890E-04, 1.300E-04, 6.470E-05, 3.600E-05, 2.250E-05 ];

Temperature = [...
  299.70,    293.70,    287.70,    283.70,    277.00, ...
  270.30,    263.60,    257.00,    250.30,    243.60, ...
  237.00,    230.10,    223.60,    217.00,    210.30, ...
  203.70,    197.00,    194.80,    198.80,    202.70, ...
  206.70,    210.70,    214.60,    217.00,    219.20, ...
  221.40,    227.00,    232.30,    237.70,    243.10, ...
  248.50,    254.00,    259.40,    264.80,    269.60, ...
  270.20,    263.40,    253.10,    236.00,    218.90, ...
  201.80,    184.80,    177.10,    177.00,    184.30, ...
  190.70,    212.00,    241.60,    299.70,    380.00];

%% ppmv
H2O = [...
  2.593E+04, 1.949E+04, 1.534E+04, 8.600E+03, 4.441E+03, ...
  3.346E+03, 2.101E+03, 1.289E+03, 7.637E+02, 4.098E+02, ...
  1.912E+02, 7.306E+01, 2.905E+01, 9.900E+00, 6.220E+00, ...
  4.000E+00, 3.000E+00, 2.900E+00, 2.750E+00, 2.600E+00, ...
  2.600E+00, 2.650E+00, 2.800E+00, 2.900E+00, 3.200E+00, ...
  3.250E+00, 3.600E+00, 4.000E+00, 4.300E+00, 4.600E+00, ...
  4.900E+00, 5.200E+00, 5.500E+00, 5.700E+00, 5.900E+00, ...
  6.000E+00, 6.000E+00, 6.000E+00, 5.400E+00, 4.500E+00, ...
  3.300E+00, 2.100E+00, 1.300E+00, 8.500E-01, 5.400E-01, ...
  4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01];
	  
CO2 = [...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,...
  3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,...
  1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01];

O3 = [...
  2.869E-02, 3.150E-02, 3.342E-02, 3.504E-02, 3.561E-02, ...
  3.767E-02, 3.989E-02, 4.223E-02, 4.471E-02, 5.000E-02, ...
  5.595E-02, 6.613E-02, 7.815E-02, 9.289E-02, 1.050E-01, ...
  1.256E-01, 1.444E-01, 2.500E-01, 5.000E-01, 9.500E-01, ...
  1.400E+00, 1.800E+00, 2.400E+00, 3.400E+00, 4.300E+00, ...
  5.400E+00, 7.800E+00, 9.300E+00, 9.850E+00, 9.700E+00, ...
  8.800E+00, 7.500E+00, 5.900E+00, 4.500E+00, 3.450E+00, ...
  2.800E+00, 1.800E+00, 1.100E+00, 6.500E-01, 3.000E-01, ...
  1.800E-01, 3.300E-01, 5.000E-01, 5.200E-01, 5.000E-01, ...
  4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04];

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


