robs      = [];
rclrSARTA = [];
rclrPCRTM = [];
rcldSARTA = [];
rcldPCRTM = [];
lat       = [];
lon       = [];
landfrac  = [];
stemp     = [];

iceamtP   = [];
wateramtP = [];
iceamtS   = [];
wateramtS = [];

icetopP   = [];
watertopP = [];
icetopS   = [];
watertopS = [];

iceszeP   = [];
waterszeP = [];
iceszeS   = [];
waterszeS = [];

dd = 01; hh = 00;
fname = ['/asl/data/rtprod_airs/2012/05/' num2str(dd,'%02d') '/xcld_ecm_41ch.airs_ctr.2012.05.'];
fname = [fname  num2str(dd,'%02d') '.' num2str(hh,'%02d') '.pcrtm.ncol50.rtp'];
      [h,ha,p,pa] = rtpread(fname);

iboo = input('Enter AIRS channel center freq : ');
dada = abs(iboo - h.vchan);
iboo = find(dada == min(dada));
fprintf(1,'closest chanID = %4i centerfreq = %8.6f \n',h.ichan(iboo),h.vchan(iboo));
oo = iboo;

iDayNight = input('enter (1) day (0) both (-1) night : ');
iLandOcean = input('enter -2 for ocean, -1 for land, 0 for land/ocean, 1:22 for TRANSCOM : ');

disp('reading in tons of files for May 2012 ...')

for dd = 1 : 31
  if mod(dd,10) == 0
    fprintf(1,'dd = %2i \n',dd)
  end
  for hh = 0 : 24
    fname = ['/asl/data/rtprod_airs/2012/05/' num2str(dd,'%02d') '/xcld_ecm_41ch.airs_ctr.2012.05.'];
    fname = [fname  num2str(dd,'%02d') '.' num2str(hh,'%02d') '.pcrtm.ncol50.rtp'];
    ee = exist(fname);
    if ee > 0
      [h,ha,p,pa] = rtpread(fname);

      find_ix

      if length(ix) > 0
        robs      = [robs      p.robs1(oo,ix)];
        rclrSARTA = [rclrSARTA p.sarta_clear(oo,ix)];
        rcldSARTA = [rcldSARTA p.rcalc(oo,ix)];
        rclrPCRTM = [rclrPCRTM p.rad_clrsky(oo,ix)];
        rcldPCRTM = [rcldPCRTM p.rad_allsky(oo,ix)];
        lat       = [lat p.rlat(ix)];
        lon       = [lon p.rlon(ix)];
        landfrac  = [landfrac p.landfrac(ix)];
        stemp     = [stemp p.stemp(ix)];

        cloud_stuff

      end
    end
  end
end      

addpath /home/sergio/MATLABCODE/SHOWSTATS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tobs = rad2bt(h.vchan(oo),robs);
tSARTA = rad2bt(h.vchan(oo),rcldSARTA);
tPCRTM = rad2bt(h.vchan(oo),rcldPCRTM);
dt = -100 : 1 : 100;
nsarta = hist(tobs - tSARTA,dt);
npcrtm = hist(tobs - tPCRTM,dt);
figure(1); plot(dt,nsarta,dt,npcrtm,'r','linewidth',2); 
  hl = legend('obs-sarta','obs-pcrtm'); set(hl,'fontsize',10); grid
  title(['all points ' num2str(h.vchan(oo)) ' cm-1']);

dt = 180 : 1 : 320;
nxobs   = hist(tobs,dt);
nxsarta = hist(tSARTA,dt);
nxpcrtm = hist(tPCRTM,dt);
figure(2); plot(dt,nxsarta,dt,nxpcrtm,'r',dt,nxobs,'k','linewidth',2); grid 
  hl = legend('sarta','pcrtm','obs'); set(hl,'fontsize',10);
  title(['all points ' num2str(h.vchan(oo)) ' cm-1']);


iProceed = input('plot histograms of comparisons between SARTA/PCRTM cloud params (-1/+1)? ');
if iProceed < 0
  error('stop')
end

figure(3); simplemap(lat,lon,tobs)
figure(4); simplemap(lat,lon,stemp-tobs)

figure(5); clf; plot(wateramtS,wateramtP,'r.',iceamtS,iceamtP,'b.'); 
  title('AMT b ice  r water'); xlabel('SARTA'); ylabel('PCRTM');
figure(6); clf; plot(waterszeS,waterszeP,'r.',iceszeS,iceszeP,'b.'); 
  title('SZE b ice  r water'); xlabel('SARTA'); ylabel('PCRTM');
figure(7); clf; plot(watertopS,watertopP,'r.',icetopS,icetopP,'b.'); 
  title('TOP b ice  r water'); xlabel('SARTA'); ylabel('PCRTM');

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); close;
figure(5)
[n x y]=myhist2d([wateramtS],[wateramtP],[[0:1:200]],[0:1:200],-1);    
subplot(121); pcolor(x,y,n); shading flat; title('AMT water'); 
xlabel('SARTA','Fontsize',12); ylabel('PCRTM','Fontsize',12); 
gaga = colormap; gaga(1,:) = 1; colormap(gaga); caxis([0 20]); colorbar('horiz');

[n x y]=myhist2d([iceamtS],[iceamtP],[[0:1:100]],[0:0.5:10],-1);    
subplot(122); pcolor(x,y,n); shading flat; title('AMT ice'); 
xlabel('SARTA','Fontsize',12); ylabel('PCRTM','Fontsize',12); 
gaga = colormap; gaga(1,:) = 1; colormap(gaga); caxis([0 20]); colorbar('horiz');

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6); close
figure(6);
[n x y]=myhist2d([waterszeS],[waterszeP],[0:1:40],[0:1:40],1);    
subplot(121); pcolor(x,y,n); shading flat; title('SZE water'); 
xlabel('SARTA','Fontsize',12); ylabel('PCRTM','Fontsize',12); 
gaga = colormap; gaga(1,:) = 1; colormap(gaga); caxis([0 20]); colorbar('horiz');

[n x y]=myhist2d([iceszeS],[iceszeP],[0:1:200],[40:0.5:80],1);    
subplot(122); pcolor(x,y,n); shading flat; title('SZE ice');
xlabel('SARTA','Fontsize',12); ylabel('PCRTM','Fontsize',12); 
gaga = colormap; gaga(1,:) = 1; colormap(gaga); caxis([0 20]); colorbar('horiz');

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7); close                                                                              
figure(7);
[n x y]=myhist2d([watertopS],[watertopP],[400:10:1000],[400:10:1000],-1);    
subplot(121); pcolor(x,y,n); shading flat; title('TOP water');
gaga = colormap; gaga(1,:) = 1; colormap(gaga); caxis([0 50]); colorbar('horiz');
xlabel('SARTA','Fontsize',12); ylabel('PCRTM','Fontsize',12); 

[n x y]=myhist2d([icetopS],[icetopP],[10:10:500],[10:10:500],-1);    
subplot(122); pcolor(x,y,n); shading flat; title('TOP ice');
gaga = colormap; gaga(1,:) = 1; colormap(gaga); caxis([0 50]); colorbar('horiz');
xlabel('SARTA','Fontsize',12); ylabel('PCRTM','Fontsize',12); 

fprintf(1,'number of points = %12i \n',length(lat))

error('lll')

%%%%%%%%%%%%%%%%%%%%%%%%%
plot(rad2bt(1231,rclrSARTA),rad2bt(1231,rclrPCRTM),'b.',rad2bt(1231,rclrSARTA),rad2bt(1231,rclrSARTA),'r.')
plot(rad2bt(1231,rclrSARTA),rad2bt(1231,rclrPCRTM) - rad2bt(1231,rclrSARTA),'.')

%[n x y]=myhist2d(fx,fy,limitsx, limitsy,LinearOrLog)
[n x y]=myhist2d(rad2bt(1231,rclrSARTA),rad2bt(1231,rclrPCRTM) - rad2bt(1231,rclrSARTA),200:1:320,-0.5:0.05:+0.5,-1);
pcolor(x,y,n); shading flat; colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%
plot(rad2bt(1231,rcldSARTA),rad2bt(1231,rcldPCRTM),'b.',rad2bt(1231,rcldSARTA),rad2bt(1231,rcldSARTA),'r.')
plot(rad2bt(1231,rcldSARTA),rad2bt(1231,rcldPCRTM) - rad2bt(1231,rcldSARTA),'.')

%[n x y]=myhist2d(fx,fy,limitsx, limitsy,LinearOrLog)
[n x y]=myhist2d(rad2bt(1231,rcldSARTA),rad2bt(1231,rcldPCRTM) - rad2bt(1231,rcldSARTA),200:1:320,-0.5:0.05:+0.5,-1);
pcolor(x,y,n); shading flat; colorbar