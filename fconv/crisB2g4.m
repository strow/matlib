
% CrIS specs for band 2 with "g4" guard channels
% 1205:0.1.250:1755 cm^-1 (441 channels)

band = 2;

v1 = 1180.0000;  % band low end {cm^-1}
v2 = 1779.9975;  % band high end {cm^-1}

vlaser = 12920;  % laser frequency {cm^-1}
% Note: require vlaser/dvc = integer
vsf = 2;         % vlaser scaling factor
rolloff = 20;

L1 = 0.4;        % longest path {cm}
Lcut = L1;       % path length saved {cm}
dvc = 1/(2*L1);  % channel spacing {cm^-1}

vcmin = 1205.00; % min channel freq
vcmax = 1755.00; % max channel freq
