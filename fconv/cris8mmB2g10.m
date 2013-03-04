
% CrIS specs for 8mm OPD band 2 with "g10" guard channels
% 1203.75:0.625:1756.25 cm^-1 (865+(10*2)=885 channels)

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

vcmin = 1203.75; % min channel freq
vcmax = 1756.25; % max channel freq
