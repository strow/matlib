
% CrIS specs for band 1
%

band = 1;

v1 = 605;        % band low end
v2 = 1190;       % band high end

L1 = 0.8;        % longest path
dvc = 1/(2*L1);  % channel spacing

vlaser = 15780;  % gives integer dvc multiples
vsf = 2;	 % vlaser scaling factor

Lcut = L1;	 % path length saved

