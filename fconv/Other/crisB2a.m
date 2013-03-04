
% CrIS specs for band 2
%

% band = 2;

v1 = 1180;        % band low end
v2 = 1780-.0025;  % band high end

% calculate interferometric parameters from the values on the 
% "CrIS Instrument Alias Unfolding" sheet, with global wlaser

global WLASER
wlaser = WLASER;
% wlaser = 775.18675;
vlaser =  1e7/wlaser;
L1 = wlaser * (10560/2) * 1e-7;
dvc = 1/(2*L1);

vsf = 2;          % vlaser scaling factor
Lcut = 1.0;       % path length saved

