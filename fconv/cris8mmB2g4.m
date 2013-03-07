
% CrIS specs for 8mm OPD band 2 with "g4" guard channels
% 1207.5:0.625:1752.5 cm^-1 (865+(4*2)=873 channels)

band = 2;

v1 = 1180.0000;  % band low end {cm^-1}
v2 = 1779.9975;  % band high end {cm^-1}

vlaser = 12920;  % laser frequency {cm^-1}
% Note: require vlaser/dvc = integer
vsf = 2;         % vlaser scaling factor
rolloff = 20;

L1 = 0.8;        % longest path {cm}
Lcut = L1;       % path length saved {cm}
dvc = 1/(2*L1);  % channel spacing {cm^-1}

vcmin = 1207.50; % min channel freq
vcmax = 1752.50; % max channel freq
