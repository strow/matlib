strow = load('two_slab_test_input.mat');
load airsheights.dat
load airslevels.dat

airsheights = flipud(airsheights);
airslevels  = flipud(airslevels);
playsN = airslevels(1:end-1)-airslevels(2:end);
playsD = log(airslevels(1:end-1)./airslevels(2:end));
airslayers = playsN./playsD;

plevs = strow.profsub.plevs;
tic; hgt1a = interp1(airslayers,airsheights,plevs,'linear'); toc
tic; hgt1b = interp1_loop(airslayers,airsheights,plevs,'linear'); toc
  sum(sum(hgt1a-hgt1b))

tic; hgt2a = interp1qr(airslayers,airsheights*ones(1,10000),plevs); toc
  sum(sum(hgt1a-hgt2a))
  imagesc((hgt1b-hgt2a)./hgt1a * 100); colormap(usa2); colorbar; caxis([-100 +100])

error('kjfa')
%% hgt2b = interp1qr2(airslayers,airsheights*ones(1,10000),plevs);
