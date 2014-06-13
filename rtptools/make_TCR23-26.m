

% Load Transcom Matrix
% 1x1 degree grid, starting at lon=-180, lat=90;
load /home/imbiriba/git/matlib/rtptools/TranscomRegionMatrix.mat

% Define new regions:

[lat lon] = mkgrid([89.5:-1:-89.5],[.5:359.5]);

RegionMatrixE = RegionMatrix;

% Add Arctic Lands as region 23
RegionMatrixE(RegionMatrixE==0 & lat> 57) = 23;
% Add Antartica as region 24
RegionMatrixE(RegionMatrixE==0 & lat<-60) = 24;
% Add Mediterranean region as region 25
RegionMatrixE(RegionMatrixE==0 & (lon>(-20+360) | lon<70) & lat>0 & lat<60) = 25;
% Add Great Lakes as region 26
RegionMatrixE(RegionMatrixE==0 & lon>(-150+360) & lon<(-50+360)) = 26;


% Save the matrix.
!mv /home/imbiriba/git/matlib/rtptools/TranscomRegionMatrix.mat /home/imbiriba/git/matlib/rtptools/TranscomRegionMatrix_orig.mat

RegionMatrix = RegionMatrixE;

save('/home/imbiriba/git/matlib/rtptools/TranscomRegionMatrix_new.mat','RegionMatrix');



