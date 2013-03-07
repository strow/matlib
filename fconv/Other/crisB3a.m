
% CrIS specs for band 3
%

% band = 3;

v1 = 2105;        % band low end
v2 = 2605-.0025;  % band high end

% calculate interferometric parameters from the values on the 
% "CrIS Instrument Alias Unfolding" sheet, with global wlaser

global WLASER
wlaser = WLASER;
% wlaser = 775.18675;
vlaser =  1e7/wlaser;
L1 = wlaser * (5200/2) * 1e-7;
dvc = 1/(2*L1);

vsf = 2;          % vlaser scaling factor
Lcut = 1.0;       % path length saved

