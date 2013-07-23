strYear = num2str(xx(1));
iMonth = xx(2);
if iMonth < 10
  strMonth = ['0' num2str(iMonth)];
else
  strMonth = [    num2str(iMonth)];
  end
iDay = xx(3);
if iDay < 10
  strDay = ['0' num2str(iDay)];
else
  strDay = [    num2str(iDay)];
  end
iGran = xx(4);
if iGran < 10
  strGran = ['00' num2str(iGran)];
elseif iGran < 100
  strGran = ['0' num2str(iGran)];
else
  strGran = [    num2str(iGran)];
  end
 
gnd = 91;
ii = find(profX.plon >= -15 & profX.plon < 45 & ...
          profX.plat >= -10 & profX.plat < 50);
plot(profX.plon(ii),profX.plat(ii),'.')
profX.plat = profX.plat(ii);
profX.plon = profX.plon(ii);
profX.T2 = profX.T2(ii);
profX.U10 = profX.U10(ii);
profX.V10 = profX.V10(ii);
profX.stemp = profX.stemp(ii);
profX.spres = profX.spres(ii);
profX.nlevs = profX.nlevs(ii);
profX.cfrac = profX.cfrac(ii);
profX.plevs = profX.plevs(:,ii);
profX.ptemp = profX.ptemp(:,ii);
profX.cc = profX.cc(:,ii);
profX.uspd = profX.uspd(:,ii);
profX.vspd = profX.vspd(:,ii);
profX.omegaspd = profX.omegaspd(:,ii);
profX.wspd = profX.wspd(:,ii);

addpath /asl/matlab/science    % usgs_deg10_dem 
[salti, landfrac] = usgs_deg10_dem(profX.plat,profX.plon); 
profX.landfrac = landfrac;

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
  newLF    = ones(1,length(newplat))*NaN;  profX.landfrac = [profX.landfrac newLF];
  newT2 = ones(1,length(newplat))*NaN;     profX.T2       = [profX.T2 newT2];
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
  newwspd  = ones(91,length(newplat))*0;   profX.wspd     = [profX.wspd newwspd];
  newospd  = ones(91,length(newplat))*0;   profX.omegaspd = [profX.omegaspd newospd];
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

figure(3); 
scatter_coast(profX.plon(ii),profX.plat(ii),150,profX.T2(ii)-profX.stemp(ii));
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
W0  = reshape(profX.wspd(91,ii),length(lon_unique),length(lat_unique));
figure(3);
ind = 1:1:240;ind = 1:5:240;
quiver(XX(ind,ind),YY(ind,ind),U0(ind,ind),V0(ind,ind)); title('gnd spd');
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
axis([-15 45 -10 50]); colorbar
figure(4)
XXX = XX(ind,ind);   XXX = XXX(:);
YYY = YY(ind,ind);   YYY = YYY(:);
ZZZ = W0(ind,ind);   ZZZ = ZZZ(:);
scatter_coast(XXX,YYY,75,ZZZ); title('gnd vertical velocity');  axis([-15 45 -10 50]); 
pause

U3km  = reshape(profX.uspd(71,ii),length(lon_unique),length(lat_unique));
V3km  = reshape(profX.vspd(71,ii),length(lon_unique),length(lat_unique));
W3km  = reshape(profX.wspd(71,ii),length(lon_unique),length(lat_unique));
figure(3);
ind = 1:1:240;ind = 1:5:240;
quiver(XX(ind,ind),YY(ind,ind),U3km(ind,ind),V3km(ind,ind)); title('2.5 km spd');
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
axis([-15 45 -10 50]); colorbar
figure(4)
ZZZ = W3km(ind,ind);   ZZZ = ZZZ(:);
scatter_coast(XXX,YYY,75,ZZZ); title('2.5 km vertical velocity'); axis([-15 45 -10 50]); 
pause

U10m  = reshape(profX.U10(ii),length(lon_unique),length(lat_unique));
V10m  = reshape(profX.V10(ii),length(lon_unique),length(lat_unique));
S10m  = U10m.^2 + V10m.^2; S10m = sqrt(S10m); 
figure(6); clf
ind = 1:1:240;ind = 1:5:240;
XXX = XX(ind,ind);   XXX = XXX(:);
YYY = YY(ind,ind);   YYY = YYY(:);
ZZZ = S10m(ind,ind); ZZZ = ZZZ(:);
hold on; scatter(XXX,YYY,100,ZZZ,'filled'); colorbar
quiver(XX(ind,ind),YY(ind,ind),U10m(ind,ind),V10m(ind,ind),'k'); title('10m spd');
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
axis([-15 45 -10 50]); 

%figure(5); contour(XX,YY,GRAD,5); colorbar; title('grad (ST)')
%hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off 

%figure(5); contour(XX,YY,SP,5); colorbar; title('spres')
%hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off 

whos XX YY GRAD ST U0 V0 U3km V3km U10m V10m T2m

saver = ...
   ['save winds_' strYear '_' strMonth '_' strDay '_' strGran ' XX YY GRAD ST U0 V0 U3km V3km U10m V10m T2m']
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
