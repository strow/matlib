
% IASI specs for 1-band spanning interferograms
%

band = 1;

v1 = 605;        % band low end {cm^-1}
v2 = 2829.9975;  % band high end {cm^-1}

L1 = 2.0;        % longest path {cm}

vlaser = 12992;  % laser frequency {cm^-1}
vsf = 2;	 % vlaser scaling factor
dvc = 0.1/(2*L1);  % channel spacing {cm^-1} x10 finer point spacing

Lcut = L1;	 % path length saved {cm}


