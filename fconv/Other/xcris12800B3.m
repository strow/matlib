
% CrIS specs for band 3
% 2155:2.500:2550 cm^-1 (159 channels)

band = 3;

v1 = 2105.0000;  % band low end {cm^-1}
v2 = 2604.9975;  % band high end {cm^-1}

vlaser = 12800;  % laser frequency {cm^-1}
% Note: require vlaser/dvc = integer
vsf = 2;         % vlaser scaling factor

L1 = 0.2;        % longest path {cm}
Lcut = L1;       % path length saved {cm}
dvc = 1/(2*L1);  % channel spacing {cm^-1}
