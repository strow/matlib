
% CrIS specs for band 1 with "g4" guard channels
% 647.5:0.625:1097.5 cm^-1 (721 channels)

band = 1;

v1 = 605.0000;   % band low end {cm^-1}
v2 = 1129.9975;  % band high end {cm^-1}

vlaser = 12920;  % laser frequency {cm^-1}
% Note: require vlaser/dvc = integer
vsf = 2;         % vlaser scaling factor
rolloff = 20;

L1 = 0.8;        % longest path {cm}
Lcut = L1;       % path length saved {cm}
dvc = 1/(2*L1);  % channel spacing {cm^-1}

vcmin =  647.5;  % min channel freq
vcmax = 1097.5;  % max channel freq
