
% Generic 1cm OPD in 15um region
% 650:0.5:2685 cm^-1 (4071 channels)

band = 1;

v1 =  630.0000;  % band low end {cm^-1}
v2 = 2704.9975;  % band high end {cm^-1}

vlaser = 12992;  % laser frequency {cm^-1}
vsf = 2;         % vlaser scaling factor

L1 = 1.0;        % longest path {cm}
Lcut = L1;       % path length saved {cm}
dvc = 1/(2*L1);  % channel spacing {cm^-1}

vcmin =  650.0;  % min channel freq
vcmax = 2685.0;  % max channel freq
