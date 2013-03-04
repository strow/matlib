
% CrIS specs for band 3
%

band = 3;

v1 = 2100;    % band low end
v2 = 2830;    % band high end

L1 = 0.2;        % longest path
dvc = 1/(2*L1);  % channel spacing

vlaser = 15780;  % gives integer dvc multiples
vsf = 2;	 % vlaser scaling factor

Lcut = 0.8;	 % path length saved

