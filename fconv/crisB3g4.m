
% CrIS specs for band 3 with "g4" guard channels
% 2145:2.500:2560 cm^-1 (167 channels)

band = 3;

v1 = 2105.0000;  % band low end {cm^-1}
v2 = 2604.9975;  % band high end {cm^-1}

vlaser = 12920;  % laser frequency {cm^-1}
% Note: require vlaser/dvc = integer
vsf = 2;         % vlaser scaling factor
rolloff = 20;

L1 = 0.2;        % longest path {cm}
Lcut = L1;       % path length saved {cm}
dvc = 1/(2*L1);  % channel spacing {cm^-1}

vcmin = 2145.00; % min channel freq
vcmax = 2560.00; % max channel freq
