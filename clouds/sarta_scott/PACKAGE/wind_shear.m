xx = input('Enter [YY MM DD GG] : ');

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

loader = ['load t2m_winds_' strYear '_' strMonth '_' strDay '_' strGran];
eval(loader);
XX = reshape(p2m.rlon,90,135);
YY = reshape(p2m.rlat,90,135);
U0  = reshape(p2m.U10,90,135);
V0  = reshape(p2m.V10,90,135);
spd10 = p2m.U10.^2 + p2m.V10.^2; spd10 = sqrt(spd10);
figure(1); clf
ctmp = coast;
hold on; scatter_coast(p2m.rlon,p2m.rlat,20,spd10);
ind1 = 1 : 5 : 90; ind2 = 1 : 5 : 135;
hold on; quiver(XX(ind1,ind2),YY(ind1,ind2),U0(ind1,ind2),V0(ind1,ind2)); 
         title('10m spd'); hold off
hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off 
axis([min(p2m.rlon) max(p2m.rlon) min(p2m.rlat) max(p2m.rlat)]); colorbar 

loader = ...
  ['load /carrot/s1/sergio/SPECIALRTPFILES/' strYear '/' strMonth '/' strDay '/'];
loader = [loader 'resetcloudparam_optimum_SST_A_Q_' strGran '_STvers6_N3_DustOnly.mat'];
eval(loader);

ii = (what2save.atrack-1)*90 + what2save.xtrack;
iDo = -1;
if iDo > 0
  p2m.plevs = p2m.plevs(:,ii);
  p2m.ptemp = p2m.ptemp(:,ii);
  p2m.uspd  = p2m.uspd(:,ii);
  p2m.vspd  = p2m.vspd(:,ii);
  p2m.T2    = p2m.T2(ii);
  p2m.U10   = p2m.U10(ii);
  p2m.V10   = p2m.V10(ii);
  p2m.stemp = p2m.stemp(ii);
  p2m.spres = p2m.spres(ii);
  p2m.rlon  = p2m.rlon(ii);
  p2m.rlat  = p2m.rlat(ii);
  end

rtpname = ['/asl/data/rtprod/' strYear '/' strMonth '/' strDay '/'];
rtpname = [rtpname 'allfov' strGran '.rtp'];
[h,ha,p,pa] = rtpread(rtpname);
if iDo > 0
  [h,p] = subset_rtp(h,p,h.glist,1:2378,ii);
  end
klayers = '/asl/packages/klayers/Bin/klayers_airs '; 
rtpwrite('dumb.ip.rtp',h,ha,p,pa);
klayerser = ['!' klayers ' fin=dumb.ip.rtp fout=dumb.op.rtp'];
eval(klayerser)
[h,ha,p,pa] = rtpread('dumb.op.rtp');
rmer = ['!/bin/rm dumb.ip.rtp dumb.op.rtp'];
eval(rmer);

if iDo > 0
  [sum(p.rlon-p2m.rlon)       sum(p.rlat-p2m.rlat) ...
   sum(p.rlon-what2save.rlon) sum(p.rlat-what2save.rlat)]
else
  [sum(p.rlon-p2m.rlon)       sum(p.rlat-p2m.rlat)]
  end

plevs     = p.plevs;
plev_hgts = p.palts;  
%[p.plevs(yuk,1) p.plays(yuk,1) p.ptemp(yuk,1) p.palts(yuk,1)]
for pp = 1:length(p.stemp)
  nlevs = p.nlevs(pp);
  plevs(nlevs,pp)     = p.spres(pp);
  plev_hgts(nlevs,pp) = p.salti(pp);
  plevs(nlevs+1:101,pp)     = -999;
  plev_hgts(nlevs+1:101,pp) = -999;
  end

figure(2); scatter_coast(what2save.rlon,what2save.rlat,50,what2save.raTau900); 
           title('\tau_{900}');

hgt = 1500;
for iii = 1 : length(ii)
  plevsx = p2m.plevs(:,ii(iii));
  dada = find(plevsx >= h2p(hgt));
  if length(dada) >= 1
    dada = dada(1);
  else 
    dada = 91;
    end
  u(iii) = p2m.uspd(dada,ii(iii));
  v(iii) = p2m.vspd(dada,ii(iii));
  end
spd = u.^2 + v.^2; spd = sqrt(spd);
figure(3); plot(spd,what2save.raTau900,'r.'); xlabel('spd'); ylabel('\tau_{900}')
figure(4); clf; scatter_coast(p2m.rlon(ii),p2m.rlat(ii),20,spd); 
  title([num2str(hgt) 'm speed']);

%%% now we can measure windshear between the levels dspeed/dz
p2m_small.plevs = p2m.plevs(:,ii);
p2m_small.ptemp = p2m.ptemp(:,ii);
p2m_small.uspd  = p2m.uspd(:,ii);
p2m_small.vspd  = p2m.vspd(:,ii);
p2m_small.T2    = p2m.T2(ii);
p2m_small.U10   = p2m.U10(ii);
p2m_small.V10   = p2m.V10(ii);
p2m_small.H10   = ones(size(p2m_small.U10))*10;  %% 10 m above ground
p2m_small.stemp = p2m.stemp(ii);
p2m_small.spres = p2m.spres(ii);
p2m_small.rlon  = p2m.rlon(ii);
p2m_small.rlat  = p2m.rlat(ii);
p2m_small.U0    = zeros(size(p2m_small.U10));    %% zero speed at gnd
p2m_small.V0    = zeros(size(p2m_small.V10));    %% zero speed at gnd
for iii = 1 : length(ii)
  thelevs = plevs(:,ii(iii));
  thehgts = plev_hgts(:,ii(iii));
  nlevs   = p.nlevs(ii(iii));
  thelevs = thelevs(1:nlevs);
  thehgts = thehgts(1:nlevs);
  palts   = interp1(log10(thelevs),thehgts,log10(p2m_small.plevs(:,iii)));
  p2m_small.palts(:,iii) = palts;
  end
p2m_small.salti = p.salti(ii);
dz = p2m_small.palts(1:90,:) - p2m_small.palts(2:91,:);
dz(91,:) = p2m_small.palts(91,:) - p2m_small.salti;
p2m_small.dz = dz;

%%% now we can do windshear; we have salti, h10,    and 91 hgts
%%%                          we have u0,v0, u10,v10 and uN,vN  N=1..91
spd   = (p2m_small.uspd).^2 + (p2m_small.vspd).^2; spd   = sqrt(spd);
spd0  = (p2m_small.U0).^2   + (p2m_small.V0).^2;   spd0  = sqrt(spd0);
spd10 = (p2m_small.U10).^2  + (p2m_small.V10).^2;  spd10 = sqrt(spd10);

shear10 = spd10./10;
shear       = spd(1:90,:) - spd(2:91,:);
shear(91,:) = spd(91,:);
shear = shear./dz;

figure(3); plot(spd10,what2save.raTau900,'.'); xlabel('10m spd'); ylabel('\tau_{900}')

figure(4); scatter_coast(p2m.rlon(ii),p2m.rlat(ii),20,spd10); title('10 m speed');
figure(4); scatter_coast(p2m.rlon(ii),p2m.rlat(ii),20,spd10/10); title('10 m shear');
figure(4); scatter_coast(p2m.rlon(ii),p2m.rlat(ii),20,shear10);    title('10 m shear');

figure(5);
scatter_coast(p2m_small.rlon,p2m_small.rlat,20,shear(91,:)); title('10 m shear')
for ll = 91 : -1 : 80
  scatter_coast(p2m_small.rlon,p2m_small.rlat,20,shear(ll,:)); title(num2str(ll));
  pause
  end
plot(shear',1:91)
