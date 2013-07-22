gnd = 91;
ii = ...
  find(profX.plon >= -15 & profX.plon < 45 & profX.plat >= -10 & profX.plat < 50);
plot(profX.plon(ii),profX.plat(ii),'.')

lon_unique = unique(profX.plon);
lat_unique = unique(profX.plat);
plot(diff(lon_unique),'.-')
plot(diff(lat_unique),'.-')
bad = find(abs(diff(lon_unique)-gridsize) > 0.01); 

%% sometimes there is no landfrac
if length(bad) > 0
  lon_unique(bad-3:bad+3)
  newplat  = lat_unique;                   profX.plat     = [profX.plat newplat];
  newspres = ones(1,length(newplat))*NaN;  profX.spres    = [profX.spres newspres];
  %newLF    = ones(1,length(newplat))*NaN;  profX.landfrac = [profX.landfrac newLF];
  newT2 = ones(1,length(newplat))*NaN;    profX.T2      = [profX.T2 newT2];
  newU10 = ones(1,length(newplat))*NaN;    profX.U10      = [profX.U10 newU10];
  newV10 = ones(1,length(newplat))*NaN;    profX.V10      = [profX.V10 newV10];
  newstemp = ones(1,length(newplat))*NaN;  profX.stemp    = [profX.stemp newstemp];
  newnlevs = ones(1,length(newplat))*0;    profX.nlevs    = [profX.nlevs newnlevs];
  newcfrac = ones(1,length(newplat))*0;    profX.cfrac    = [profX.cfrac newcfrac];
  newplevs = ones(91,length(newplat))*0;   profX.plevs    = [profX.plevs newplevs];
  newptemp = ones(91,length(newplat))*0;   profX.ptemp    = [profX.ptemp newptemp];
  newcc    = ones(91,length(newplat))*0;   profX.cc       = [profX.cc newcc];
  newuspd  = ones(91,length(newplat))*0;   profX.uspd     = [profX.uspd newuspd];
  newvspd  = ones(91,length(newplat))*0;   profX.vspd     = [profX.vspd newvspd];
  newplon  = ones(1,length(newplat))*(lon_unique(bad)+0.5);
     profX.plon     = [profX.plon newplon];
  end

ii = ...
  find(profX.plon >= -15 & profX.plon < 45 & profX.plat >= -10 & profX.plat < 50);
plot(profX.plon(ii),profX.plat(ii),'.')
lon_unique = unique(profX.plon(ii));
lat_unique = unique(profX.plat(ii));

ctmp = coast; 
figure(1); scatter_coast(profX.plon(ii),profX.plat(ii),150,profX.stemp(ii));
title('stemp')

figure(2); scatter_coast(profX.plon(ii),profX.plat(ii),150,profX.T2(ii));
title('T2m')

figure(3); scatter_coast(profX.plon(ii),profX.plat(ii),150,profX.T2(ii)-profX.stemp(ii));
title('T2m-stemp')

pause

T2m = reshape(profX.T2(ii),length(lon_unique),length(lat_unique));
ST  = reshape(profX.stemp(ii),length(lon_unique),length(lat_unique));
SP  = reshape(profX.spres(ii),length(lon_unique),length(lat_unique));
[FX,FY] = gradient(ST);
GRAD = sqrt(FX.*FX + FY.*FY);
figure(2); scatter_coast(profX.plon(ii),profX.plat(ii),30,GRAD(:)');
title('stemp gradient');

XX = reshape(profX.plon(ii),length(lon_unique),length(lat_unique));
YY = reshape(profX.plat(ii),length(lon_unique),length(lat_unique));
U0  = reshape(profX.uspd(91,ii),length(lon_unique),length(lat_unique));
V0  = reshape(profX.vspd(91,ii),length(lon_unique),length(lat_unique));
figure(3);
ind = 1:1:120;ind = 1:5:120;
quiver(XX(ind,ind),YY(ind,ind),U0(ind,ind),V0(ind,ind)); title('gnd spd');
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
axis([-15 45 -10 50]); colorbar

U3km  = reshape(profX.uspd(71,ii),length(lon_unique),length(lat_unique));
V3km  = reshape(profX.vspd(71,ii),length(lon_unique),length(lat_unique));
figure(4);
ind = 1:1:120;ind = 1:5:120;
quiver(XX(ind,ind),YY(ind,ind),U3km(ind,ind),V3km(ind,ind)); title('2.5 km spd');
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
axis([-15 45 -10 50]); colorbar

U10m  = reshape(profX.U10(ii),length(lon_unique),length(lat_unique));
V10m  = reshape(profX.V10(ii),length(lon_unique),length(lat_unique));
figure(5);
ind = 1:1:120;ind = 1:5:120;
quiver(XX(ind,ind),YY(ind,ind),U10m(ind,ind),V10m(ind,ind)); title('10m spd');
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
axis([-15 45 -10 50]); colorbar

figure(5); contour(XX,YY,GRAD,5); colorbar; title('grad (ST)')
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off 

figure(5); contour(XX,YY,SP,5); colorbar; title('spres')
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off 

whos XX YY GRAD ST U0 V0 U3km V3km U10m V10m

iYear = num2str(xx(1));
iMonth = xx(2);
if iMonth < 10
  iMonth = ['0' num2str(iMonth)];
else
  iMonth = [    num2str(iMonth)];
  end
iDay = xx(3);
if iDay < 10
  iDay = ['0' num2str(iDay)];
else
  iDay = [    num2str(iDay)];
  end
iGran = xx(4);
if iGran < 10
  iGran = ['00' num2str(iGran)];
elseif iGran < 100
  iGran = ['0' num2str(iGran)];
else
  iGran = [    num2str(iGran)];
  end

saver = ...
   ['save winds_' iYear '_' iMonth '_' iDay '_' iGran ' XX YY GRAD ST U0 V0 U3km V3km U10m V10m T2m']
eval(saver)