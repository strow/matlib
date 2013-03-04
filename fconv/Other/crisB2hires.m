
% CrIS specs for band 2
%

band = 2;

v1 = 1200;       % band low end
v2 = 1800;       % band high end

L1 = 0.8;        % longest path
dvc = 1/(2*L1);  % channel spacing

vlaser = 15780;  % gives integer dvc multiples
vsf = 2;	 % vlaser scaling factor

Lcut = L1;	 % path length saved

