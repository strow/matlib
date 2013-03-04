
% TANSO TIR specs for 1-band spanning interferograms

band = 1;

v1 = 630;        % band low end {cm^-1}
v2 = 1879.9975;  % band high end {cm^-1}

L1 = 2.51470464; % longest path {cm}

vlaser = 12992;  % laser frequency {cm^-1}
vsf = 2;	 % vlaser scaling factor
dvc = 1/(2*L1);  % channel spacing {cm^-1}

Lcut = L1;	 % path length saved {cm}

rolloff = 20;    % Rolloff spectral distance {cm^-1}

