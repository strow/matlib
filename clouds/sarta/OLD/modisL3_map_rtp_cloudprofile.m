function waterDME_MODIS = modisL3_map_rtp_cloudprofile(prof);

%% same as modisL3_map_rtp except this is for all rlat/rlon in the profile sructure, as 
%% opposed to only where there are ctype/ctype2 = 101

addpath /asl/matlib/aslutil/

[yy,mm,dd] = tai2utc(prof.rtime);

yy = ceil(mean(yy));
mm = ceil(mean(mm));
fname = ['/asl/s1/sergio/MODIS_L3/modisL3dme_' num2str(yy) '_' num2str(mm,'%02d') '.mat'];
if ~exist(fname)
  fname = ['/asl/s1/sergio/MODIS_L3/modisL3dme_2011_03.mat'];
  disp('loading in WRONG modis L3 file modisL3dme_2011_03.mat')
end
loader = ['load('''  fname ''');'];
eval(loader);

cmap_modis_michaelking = [43 45 44; 58 49 73; 52 65 82; 34 77 93; 57 81 44; 94 85 19; 
                          99 62 19; 94 33 13; 80 07 15; 36 16 21];
cmap_modis_michaelking = cmap_modis_michaelking/100;

%figure(3); scatter(rlon,-rlat,10,cfrac(:)); title('cfrac')

figure(4); scatter(rlon(:),-rlat(:),20,waterDME(:)); title('water rme'); colorbar
  caxis(2*[5 25]);; colormap(cmap_modis_michaelking)
pause(0.1)
figure(5); scatter(rlon(:),-rlat(:),100,waterDME(:)); title('water rme'); colorbar
  caxis(2*[5 25]);; colormap(cmap_modis_michaelking)
pause(0.1)
axis([min(prof.rlon)-2 max(prof.rlon)+2 min(prof.rlat)-2 max(prof.rlat)+2])

%figure(5); scatter(rlon,-rlat,10,iceDME(:)); title('ice rme')
%  caxis([16 32]);; colormap(cmap_modis_michaelking)

ii = sub2ind(size(waterDME),max(1,min(180,round(-prof.rlat+89.5))),mod(round(prof.rlon-0.5)+180,360)+1);
waterDME_MODIS = waterDME(ii);
figure(1); scatter(prof.rlon,prof.rlat,30,waterDME_MODIS)
caxis(2*[5 25]);; colormap(cmap_modis_michaelking); title('MODIS dme'); colorbar

