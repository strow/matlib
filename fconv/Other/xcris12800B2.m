
% CrIS specs for band 2
% 1210:0.1.250:1750 cm^-1 (433 channels)

band = 2;

v1 = 1180.0000;  % band low end {cm^-1}
v2 = 1779.9975;  % band high end {cm^-1}

vlaser = 12800;  % laser frequency {cm^-1}
% Note: require vlaser/dvc = integer
vsf = 2;         % vlaser scaling factor

L1 = 0.4;        % longest path {cm}
Lcut = L1;       % path length saved {cm}
dvc = 1/(2*L1);  % channel spacing {cm^-1}
