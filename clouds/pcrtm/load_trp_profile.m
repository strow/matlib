%% see AFGL proflies in
%% DFAFGL in /asl/packages/klayersV205/Src/incLAY.f is
%%    /asl/packages/klayersV205/Data/glatm.dat -> glatm_16Aug2010.dat
%% [h,ha,p,pa] = rtpread('/asl/packages/klayersV205/Data/adafgl_16Aug2010_ip.rtp');

disp('    loading in TRP profile')

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

CH4 = [...
  1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, ...
  1.700E+00, 1.700E+00, 1.699E+00, 1.697E+00, 1.693E+00, ...
  1.685E+00, 1.675E+00, 1.662E+00, 1.645E+00, 1.626E+00, ...
  1.605E+00, 1.582E+00, 1.553E+00, 1.521E+00, 1.480E+00, ...
  1.424E+00, 1.355E+00, 1.272E+00, 1.191E+00, 1.118E+00, ...
  1.055E+00, 9.870E-01, 9.136E-01, 8.300E-01, 7.460E-01, ...
  6.618E-01, 5.638E-01, 4.614E-01, 3.631E-01, 2.773E-01, ...
  2.100E-01, 1.651E-01, 1.500E-01, 1.500E-01, 1.500E-01, ...
  1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01, ...
  1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02];

