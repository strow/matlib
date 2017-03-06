%{

http://www.ecmwf.int/sites/default/files/elibrary/2011/9441-new-prognostic-bulk-microphysics-scheme-ifs.pdf

A new prognostic bulk microphysics scheme for the IFS Richard
M. Forbes 1 , Adrian M. Tompkins 2 and Agathe Untch 1 Research
Department 1 ECMWF 2 ICTP, Italy September 2011
Report 649
Page 12

In the IFS model the ice particle effective radius is a function of
temperature and ice water content, based on analysis of observational
aircraft data from Sun and Rikus ( 1999 ) (revised by Sun , 2001 ),
which covers a range of particle sizes.  The effective radius is
currently limited between a minimum of 20+40cos( latitude ) μ m and a
maximum of 155 μ m.  A factor of 0.64952 is used to convert from ef-
fective radius to particle diameter.  At present, the same op tical
properties are used for the ice and snow particles; Fu ( 1996 ) for
the shortwave optical properties and Fu et al.  ( 1998 ) for the
longwave spectral emissivity.

Sun, Z. (2001). Reply to comments by Greg M. McFarquhar on ”Pa
rametrization of effective sizes of cirrus-cloud particles and its
verification against observ ations”.  Q. J. R. Meteorol. Soc.  , 127 ,
267–271.

Sun, Z. and Rikus, L. (1999). Parametrization of effective s izes of
cirrus-cloud particles and its verifica- tion against observations.
Q. J. R. Meteorol. Soc.  , 125 , 3037–3055

Eqns 9,10,11

a linear fit of the effective size, De,
versus temperature for a given ice-water content.
The equation has the form,

De = a(IWC) + b(IWC)(T + 190), (9)

where T is the temperature in Celsius and De is in um.  The
coefficients a(IWC) and b(IWC) are further fitted to IWC using the
following expressions,

a(IWC) = 26.1571 |1ogl0(10^-12 + IWC)|^(-0.5995)   (10)
b(IWC) = 0.6402 + 0.1810 log10(10^-12 + IWC)       (11) 

%%%%%%%%%%%%%%%%%%%%%%%%%

First ice cloud effective particle size parameterization based on
combined lidar and radar data D. P. Donovan and A. C. Lamneren, GRL
VOL.  29, NO.  1, 1006, 10.1029/2001GL013731, 2002

r eff = sum(i=0,2) Ai T^i
                      A0    A1     A2
Complex Polycrystals 356.9  2.65  -0.0374
Bullet Rossetts      718.6  5.66  -0.0705

T = 200:273;
T = T-273;
cp = 356.9 + 2.65*T -0.0374*T.^2;
br = 718.6 + 5.66*T -0.0705*T.^2;
plot(T+273,cp,T+273,br)

%}


%% see eqn 2.18, pg 11 of  ECMWF   9211-part-iv-physical-processes.pdf
%% looks like 1 = TOA
%%            N = GND

addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/SHOWSTATS
addpath /home/sergio/MATLABCODE/LOADMIE
addpath /home/sergio/MATLABCODE/PLOTTER

if ~exist('p')
  [h,ha,p,pa] = rtpread('//asl/rtp/rtprod_airs/2011/03/11/ECM_CLOUD_SARTA3/cloudy_airs_l1b_ecm_sarta3_baum_ice_modiswater.2011.03.11.039.rtp');

  %figure(1); plot(p.sarta_lvlODice(:),p.ciwc(:),'.'); xlabel('ice OD'); ylabel('ciwc')
  %figure(2); plot(p.sarta_lvlODwater(:),p.clwc(:),'.'); xlabel('water OD'); ylabel('clwc')

  jett = jet; jett(1,:) = 1;
  figure(1); [n,nx,ny] = myhist2d(p.sarta_lvlODice(:),p.ciwc(:),0:0.025:10,0:1e-5:5e-4);
    pcolor(0:0.025:10,0:1e-5:5e-4,log10(n)); xlabel('ice OD'); ylabel('ciwc'); shading flat; colorbar; colormap(jett)

  figure(2); [n,nx,ny] = myhist2d(p.sarta_lvlODwater(:),p.clwc(:),0:0.025:10,0:1e-6:1e-4);
    pcolor(0:0.025:10,0:1e-6:1e-4,log10(n)); xlabel('water OD'); ylabel('clwc'); shading flat; colorbar; colormap(jett)  
  figure(2); [n,nx,ny] = myhist2d(sum(p.sarta_lvlODwater,1),sum(p.clwc,1),0:5:500,0:1e-4:1.5e-3);  
    pcolor(0:5:500,0:1e-4:1.5e-3,log10(n)); xlabel('water OD'); ylabel('clwc'); shading flat; colorbar; colormap(jett)

  %% for ice   looks like OD 5   corresponds to 15e-5 g/g
  %% for water looks like OD 200 corresponds to 0.5e-3 g/g
end

if ~exist('p0')
  load /asl/rtp/rtprod_airs/2009/03/01/forITOVS_ECM.mat
end

if ~exist('p1ALL') & ~exist('p1ALL_7377_2009_03_01.mat')
  addpath /home/sergio/MATLABCODE/matlib/clouds/pcrtm/
  addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
  run_sarta.clear = -999;
  run_sarta.cloud = -999;
  [h1ALL,h1aALL,p1ALL,p1aALL] = driver_pcrtm_cloud_rtp(h,ha,p0,pa,run_sarta);
  p1ALL = rmfield(p1ALL,'rcalc');
  p1ALL = rmfield(p1ALL,'rad_allsky_std');
  p1ALL = rmfield(p1ALL,'rad_clrsky');
  p1ALL = rmfield(p1ALL,'rad_allsky');  
  p1ALL = rmfield(p1ALL,'rcalc_std');
  p1ALL = rmfield(p1ALL,'sarta_clear');
  p1ALL = rmfield(p1ALL,'sarta_cloud');
  p1ALL = rmfield(p1ALL,'sarta_cloud_cumsumMinus1');        
  save p1ALL_7377_2009_03_01.mat h1ALL h1aALL p1ALL p1aALL run_sarta

  %% two slab clouds
  fip = 'forITOVS_ECM_2slab.ip.rtp';
  fop = 'forITOVS_ECM_2slab.op.rtp';
  rtpwrite(fip,h1ALL,h1aALL,p1ALL,p1aALL);  
  klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
  klayerser = ['!' klayers ' fin=' fip ' fout=' fop];
  eval(klayerser)
  
  %% see /home/sergio/MATLABCODE/matlib/clouds/sarta/get_sarta_cloud100layer.m
  fip = 'forITOVS_ECM_100layercloud.ip.rtp';
  fop = 'forITOVS_ECM_100layercloud.op.rtp';
  gas_str = 'nwant=10 listg=1,2,3,4,5,6,9,12,201,202 ';  
  h1ALL_ip = h1ALL;
  h1ALL_ip.ngas = 8;
  h1ALL_ip.glist = [h1ALL.glist; [201 202]'];
  h1ALL_ip.gunit = [h1ALL.gunit; [21  21]'];
  p1ALL_ip = p1ALL;
  p1ALL_ip.gas_201 = p1ALL.clwc;  %% 201 = water cloud ~ ctype 101
  p1ALL_ip.gas_202 = p1ALL.ciwc;  %% 202 = ice cloud   ~ ctype 201
  rtpwrite(fip,h1ALL_ip,h1aALL,p1ALL_ip,p1aALL);  
  klayers100 = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';
  klayerser = ['!' klayers100 ' fin=' fip ' fout=' fop ' ' gas_str];
  eval(klayerser)
  [h1ALL_op,h1aALL,p1ALL_op,p1aALL] = rtpread(fop);
  
elseif ~exist('p1ALL') & exist('p1ALL_7377_2009_03_01.mat')
  load p1ALL_7377_2009_03_01.mat
  fop = 'forITOVS_ECM_100layercloud.op.rtp';
  [h1ALL,h1aALL,p1ALL_op,p1aALL] = rtpread(fop);  
end


%% pcrtm_ODx is ALL layers, pcrtm_OD is only lower layers ... so latter is larger than former
x = p1ALL.pcrtm_iceODX; y = p1ALL.pcrtm_iceOD; plot(x,mean(p1ALL.pcrtm_lvlODice),'.',y,mean(p1ALL.pcrtm_lvlODice),'.')
dn = 0:1e3:1e6; semilogy(dn,hist(p1ALL.pcrtm_lvlODice(:) ./ (eps+p1ALL.ciwc(:)),dn))

%% from Sun/Rikus 1991 QJRMS
%% need IWC in g/m3 : see /home/sergio/MATLABCODE/matlib/clouds/pcrtm/PCRTM_compute_for_AIRS_spectra.m line 310
R = 287.05;

%% test
qi = -4 : 0.1 : 0; qi = 10.^qi;
a = 26.1571 * abs(log10(10^(-12) + qi)).^(-0.5995);
b = 0.6402  + 0.1810 * (log10(10^(-12) + qi));
test1 =  a + b*(-10 + 190);
test2 =  a + b*(-60 + 190);
plot(qi,test1,qi,test2,'linewidth',2); set(gca,'xscale','log'); axis([1e-5 1 0 150]); grid
hl = legend('-10C','-60C','location','northwest');
set(gca,'YTick',[0:10:150])

qi = p1ALL.ciwc ./ p1ALL.ptemp  .* p1ALL.plevs * 100/R *1e3;
qi = max(eps,qi);
qw = p1ALL.clwc ./ p1ALL.ptemp  .* p1ALL.plevs * 100/R *1e3;
qw = max(eps,qw);
a = 26.1571 * abs(log10(10^(-12) + qi)).^(-0.5995);
b = 0.6402  + 0.1810 * (log10(10^(-12) + qi));
ice_deff_sun_rikus = a + b.*(p1ALL.ptemp-273 + 190);;
woo = find(p1ALL.ciwc < eps | imag(ice_deff_sun_rikus) > eps | real(ice_deff_sun_rikus) < eps); ice_deff_sun_rikus(woo) = NaN;

%%% >>>>>>>>>>>>>>>>>>>>>>>>>
playsN = p1ALL_op.plevs(1:100,:)-p1ALL_op.plevs(2:101,:);
playsD = log(p1ALL_op.plevs(1:100,:)./p1ALL_op.plevs(2:101,:));
p1ALL_op.plays = [playsN./playsD; zeros(1,7377)];
plot(p1ALL_op.gas_201(:,1),p1ALL_op.plays(:,1),'b',qw(:,1)*287,p1ALL.plevs(:,1),'c',...
     p1ALL_op.gas_202(:,1),p1ALL_op.plays(:,1),'r',qi(:,1)*287,p1ALL.plevs(:,1),'m','linewidth',2)
%% klayers output in g/m2 while ciwc,clwc in g/m3
dz = [p1ALL_op.palts(1:100,:) - p1ALL_op.palts(2:101,:); zeros(1,7377)];
ix = 1;
ix = 100;
plot(p1ALL_op.gas_201(:,ix)./dz(:,ix),p1ALL_op.plays(:,ix),'b',qw(:,ix),p1ALL.plevs(:,ix),'c',...
     p1ALL_op.gas_202(:,ix)./dz(:,ix),p1ALL_op.plays(:,ix),'r',qi(:,ix),p1ALL.plevs(:,ix),'m','linewidth',2)
set(gca,'ydir','reverse')
%%% >>>>>>>>>>>>>>>>>>>>>>>>>

% coefficents from S-C Ou, K-N. Liou, Atmospheric Research
% 35(1995):127-138.
% for computing ice cloud effective size
c0 = 326.3;
c1 = 12.42;
c2 = 0.197;
c3 = 0.0012;
tcld = p1ALL.ptemp - 273.16;
tcld = min(max(tcld,-50),-25);
cldde_ice_ou = c0 + c1 * tcld + c2 * tcld.^2 + c3 * tcld.^3;
woo = find(p1ALL.ciwc < eps); cldde_ice_ou(woo) = NaN;

plot(p1ALL.rlat,ice_deff_sun_rikus,'b.',p1ALL.rlat,cldde_ice_ou,'r.',...
     p1ALL.rlat,2*(20+40*cos(p1ALL.rlat*pi/180)),'k.',p1ALL.rlat,ones(size(p1ALL.rlat))*155*2,'k.')
axis([-90 +90 0 200])


figure(4)
  dn = 0:200; plot(dn,hist(ice_deff_sun_rikus(:),dn),'b',dn,hist(cldde_ice_ou(:),dn),'r')
  dn = 0:200; plot(dn,hist(ice_deff_sun_rikus(:),dn),'b',dn,hist(nanmean(ice_deff_sun_rikus),dn),'c',...
                   dn,hist(cldde_ice_ou(:),dn),'r',dn,hist(p1ALL.pcrtm_iceDME,dn),'m')
  set(gca,'yscale','log'); xlabel('Deff'); ylabel('Count')
figure(5)
  plot(1:7377,nanmean(ice_deff_sun_rikus),1:7377,p1ALL.pcrtm_iceDME);
  ylabel('Deff (um)'); title('(b) sun rikus ECMWF (r) ou SARTA/PCRTM')

if ~exist('Ei')
  iceDME = 10 * ones(size(p1ALL.pcrtm_iceDME));
  yes = find(isfinite(p1ALL.pcrtm_iceDME) & p1ALL.pcrtm_iceDME > 0);
  no  = find(isnan(p1ALL.pcrtm_iceDME) | isinf(p1ALL.pcrtm_iceDME) |  p1ALL.pcrtm_iceDME < eps);
  iceDME(yes) = p1ALL.pcrtm_iceDME(yes);

  iceDME = 10 * ones(size(p1ALL.pcrtm_iceDME));
  sun_rikus_dme = nanmean(ice_deff_sun_rikus);
  yes = find(isfinite(sun_rikus_dme) & sun_rikus_dme > 0);
  no  = find(isnan(sun_rikus_dme) | isinf(sun_rikus_dme) |  sun_rikus_dme < eps);
  iceDME(yes) = sun_rikus_dme(yes);

  [Ei,Wi,Gi] = loadmie_tau_EWG(960,iceDME  ,'I',4,-1,-1,ones(size(p1ALL.pcrtm_iceDME)),1);
  Ei(no) = 0;
  Wi(no) = 0;
  Gi(no) = 0;
  ice_yes = ones(size(p1ALL.pcrtm_iceDME));
  ice_yes(no) = 0;
  
  waterOD = ones(size(p1ALL.pcrtm_iceOD)) * 20;
  [Ew,Ww,Gw] = loadmie_tau_EWG(960,waterOD,'W',250,1,6,ones(size(p1ALL.pcrtm_iceOD)),1);
  no  = find(isnan(p1ALL.pcrtm_waterODX) | isinf(p1ALL.pcrtm_waterODX) | p1ALL.pcrtm_waterODX < eps);
  Ew(no) = 0;
  Ww(no) = 0;
  Gw(no) = 0;
  water_yes = ones(size(p1ALL.pcrtm_iceOD));
  water_yes(no) = 0;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
(Andy Tangborn) I found a paper on the ECMWF cloud physics (maximum random
overlap; Jakob, et al , 1999 in QJRMS).  They use the following
algorithm (where Ck is the cloud fraction from the kth layer to the
TOA, ak is the cloud fraction for a single layer k):
  Ck = 1 -(1-C_{k-1} [ 1 - max(a_{k-1},a_{k})]/[1 - min(a_{k-1}, 1-delta)]
  for k=2:n_levels.
  init cond : for k=1, c1=a1.
  where delta = 0.000001
%}

delta = 0.000001;
for nn = 1 : p0.nlevs(1)
  an = p0.cc(nn,:);
  if nn == 1
    cn = an;
  else
    an_1 = p0.cc(nn-1,:);
    bnum = 1 - max(an,an_1);
    bden = 1 - min(an_1,1-delta);
    cn = 1 - (1-cn).*bnum./bden;
  end
end
tcc_andy = cn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (A) this DOES NOT include exponential factor for cloud overlap
for pp = 1 : length(p0.stemp)
  cc = p0.cc(:,pp);
  cc = p0.cc(2:92,pp);  
  oneminus = 1-cc;
  gah = cumprod(oneminus);
  tcc_sergio(pp) = 1 - gah(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
ccx = ones(size(cc));  %% test
ccx = zeros(size(cc));  %% test
ccx = 1-rand(size(cc));  %% test
ccx = cc;  %% from the last profile

oneminus = 1-ccx;
gah = cumprod(oneminus);
figure(1); plot(ccx,'o-'); title('cc');
figure(2); plot(oneminus,'o-'); title('1-cc')
figure(3); plot(cumprod(oneminus),'o-'); title('cumprod(oneminus)'); fprintf(1,'tcc_sergio = %8.6f \n',1-gah(end))

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); plot(p0.tcc,tcc_sergio,'.'); line([0 1],[0 1],'color','r')
  xlabel('TCC ECMWF'); ; ylabel('TCC SERGIO');
figure(5); plot(p0.tcc,p0.tcc-tcc_sergio,'.'); xlabel('TCC ECMWF'); ; ylabel('\delta TCC (ECMWF-SERGIO)');
figure(6); dn = -1 : 0.025 : +1; plot(dn,hist(p0.tcc-tcc_sergio,dn),dn,hist(p0.tcc-tcc_andy,dn),'linewidth',2)
  legend('sergio','andy'); xlabel('ECMWF - attempts')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   using rough factors %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gg_to_od = input('Enter rough conversion from g/g to OD (0 makes TCC too small, 1e5 ok, 1e60 makes TCC like above) : NOT USED!!!! : ');

clear *fac*

%% (B) this DOES include exponential factor for cloud overlap, eqn 2.19
%% which depends on effective asymmetry, SSA and optical depth
%% so it depends on cloud amount
for pp = 1 : length(p0.stemp)
  cc = p0.cc(:,pp);
  cc = p0.cc(2:92,pp);
  ciwc = p0.ciwc(2:92,pp);
  clwc = p0.clwc(2:92,pp);
  
  fac = gg_to_od*(ciwc + clwc);            %%% total OD -- 1000 is a fudge to go from g/g to OD
  fac = ciwc*5/15e-5 + clwc*200/0.5e-3;    %%% total OD -- 1000 is a fudge to go from g/g to OD  
  fac = (1-0.5*1*1)*fac;                   %% pretend SSA ~ 0.5, g ~ 1
  expfac = 1 - exp(-fac);

  oneminus = 1-cc.*expfac;
  gah = cumprod(oneminus);
  tcc(pp) = 1 - gah(end);
end

%% (C) this DOES include exponential factor for cloud overlap, eqn 2.19
%% which depends on effective asymmetry, SSA and optical depth
%% so it depends on cloud amount, use more accurate ODs

fudge = 0.05;   %% modify eqn in ECMWF Tech report
fudge = 1;      %% eqn in ECMWF Tech report

fudge

for pp = 1 : length(p0.stemp)
  cc = p0.cc(:,pp);
  cc = p0.cc(2:92,pp);
  ciwc = p0.ciwc(2:92,pp);
  clwc = p0.clwc(2:92,pp);

  %% simple try, constant SSA and asymm
  odfac  = p1ALL.pcrtm_lvlODice(2:92,pp)*Ei(pp) + p1ALL.pcrtm_lvlODwater(2:92,pp)*Ew(pp); %% total OD
  fac    = (1-0.5*1*1)*odfac;  %% pretend SSA ~ 0.5, g ~ 1
  expfac = 1 - exp(-fac);

  %% more involved try, to use better ice SSA and sym (and ditto for water)
  %% we actually need atmospheric OD, but forget this for now
  od_level_i  = p1ALL.pcrtm_lvlODice(2:92,pp) + p1ALL.pcrtm_lvlODwater(2:92,pp); %% total OD
  od_level_i  = p1ALL.pcrtm_lvlODice(2:92,pp)*ice_yes(pp) + p1ALL.pcrtm_lvlODwater(2:92,pp)*water_yes(pp); %% total OD  
  ssa_level_i = (Wi(pp) + Ww(pp))/2;
  asym_level_i = Gi(pp)*mean(p1ALL.pcrtm_lvlODice(2:92,pp)) + Gw(pp)*mean(p1ALL.pcrtm_lvlODwater(2:92,pp));
  asym_level_i = asym_level_i / (mean(p1ALL.pcrtm_lvlODice(2:92,pp)) + mean(p1ALL.pcrtm_lvlODwater(2:92,pp)));
  fac    = (1-ssa_level_i.*asym_level_i.*asym_level_i).*od_level_i;
  expfac = 1 - exp(-fac);

  oneminus = 1-cc.*expfac;
  gah = cumprod(oneminus);
  tcc_better(pp) = 1 - gah(end);
end

plot(1-cumprod(oneminus),1:91,'o-'); set(gca,'ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%
ccx = ones(size(cc));  %% test
ccx = zeros(size(cc));  %% test
ccx = 1-rand(size(cc));  %% test
ccx = cc;  %% from the last profile

honeminus = 1-ccx;
hah = cumprod(honeminus);
figure(1); plot(ccx,'o-'); title('cc');
figure(2); plot(honeminus,'o-'); title('1-cc')
figure(3); plot(cumprod(honeminus),'o-'); title('cumprod(honeminus)'); fprintf(1,'tcc = %8.6f \n',1-hah(end))
figure(4); plot(1:91,ciwc,1:91,clwc);
figure(5); plot(expfac);

%%%%%%%%%%%%%%%%%%%%%%%%%
tcc_andy(tcc_andy > 1) = 1;

figure(6); plot(p0.tcc,tcc,'m.',p0.tcc,tcc_better,'r.',p0.tcc,tcc_sergio,'b.'); line([0 1],[0 1],'color','k')
  xlabel('TCC ECMWF'); ; hl = legend('TCC ECMWF weighted CIWC,CLWC','TCC ECMWF lvl OD weighted CIWC,CLWC','TCC SERGIO simple no weight');
  set(hl,'fontsize',10);  
figure(7); plot(p0.tcc,p0.tcc-tcc,'m.',p0.tcc,p0.tcc-tcc_better,'r.',p0.tcc,p0.tcc-tcc_sergio,'b.',p0.tcc,p0.tcc-tcc_andy,'k.');
  xlabel('TCC ECMWF'); ; ylabel('\delta TCC (ECMWF-SERGIO)');
  hl = legend('weighted CIWC,CLWC','weighted CIWC,CLWC 2','simple no weight','andy 1999 QJRMS'); set(hl,'fontsize',10);
figure(8); dn = -1 : 0.025 : +1;
  lenP = length(p0.stemp);
  semilogy(dn,hist(p0.tcc-tcc,dn)/lenP,'m',dn,hist(p0.tcc-tcc_better,dn)/lenP,'r',dn,hist(p0.tcc-tcc_sergio,dn)/lenP,'b',dn,hist(p0.tcc-tcc_andy,dn)/lenP,'k','linewidth',2)
  hold on
  ocean = find(p0.landfrac == 0);
  semilogy(dn,hist(p0.tcc(ocean)-tcc(ocean),dn)/lenP,'m--',...
           dn,hist(p0.tcc(ocean)-tcc_better(ocean),dn)/lenP,'r--',...
	   dn,hist(p0.tcc(ocean)-tcc_sergio(ocean),dn)/lenP,'b--',...
	   dn,hist(p0.tcc(ocean)-tcc_andy(ocean),dn)/lenP,'k--',...	   
	   'linewidth',2)
  hold off
  hl = legend('weighted CIWC,CLWC','weighted LVL OD CIWC, CLWC','simple no weight','andy 1999 QJRMS'); set(hl,'fontsize',10);
  xlabel('TCC from ECMWF - calcs');
  title('Solid : land/ocean; dashed : ocean only')
  grid

figure(9); scatter_coast(p0.rlon,p0.rlat,10,p0.tcc);               colormap jet; title('p0.tcc')
figure(10); scatter_coast(p0.rlon,p0.rlat,10,tcc_better);          colormap jet; title('weighted LVL OD CIWC, CLWC --> tcc')
figure(11); scatter_coast(p0.rlon,p0.rlat,10,p0.tcc - tcc_better); colormap jet; title('p0.tcc - weighted LVL OD CIWC, CLWC --> tcc'); caxis([-0.25 +0.25])
figure(12); plot(p0.tcc-tcc_better,p0.landfrac,'.')
figure(12); scatter(p0.tcc,p0.tcc-tcc_better,10,p0.landfrac,'filled'); colorbar; colormap jet

figure(9)