addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[h0,ha,p0,pa] = rtpread('/asl/s1/sergio/LES/RTP/rico.50446.0.nc_perturbWV_10.rtp');
%[h0,ha,p0,pa] = rtpread('/asl/data/rtprod_airs/2011/03/11/cloudy_airs_l1b_era_sarta_baran_ice.2011.03.11.240.rtp');
%[h0,ha,p0,pa] = rtpread('/asl/data/rtprod_airs/2011/03/11/cloudy_airs_l1b_era_sarta_baum_ice.2011.03.11.240.rtp');
[h0,ha,p0,pa] = rtpread('/asl/data/rtprod_airs/2011/03/11/cloudy_airs_l1b_era_sarta_baum_ice.2011.03.11.230.rtp');
%[h0,ha,p0,pa] = rtpread('/asl/data/rtprod_airs/2011/03/11/cloudy_airs_l1b_era_sarta_baum_ice.2011.03.11.001.rtp');

[h0,ha,p0,pa] = rtpadd_emis_wis(h0,ha,p0,pa);
iJump = 1000;
iJump = 100;
iJump = 25;
iJump = 10;
[h0,p0] = subset_rtp_clouds(h0,p0,[],[],1:iJump:length(p0.stemp));

addpath /home/sergio/MATLABCODE/matlib/clouds/sarta/
addpath /home/sergio/MATLABCODE/matlib/clouds/pcrtm/
  run_sarta.clear = +1;
  run_sarta.cloud = +1;
    run_sarta.ncol  = 50;
    run_sarta.ncol  = 25;
  run_sarta.cfrac = 1;
  run_sarta.ice_water_separator = 440;

  tstart = tic;
  p1 = driver_sarta_cloud_rtp(h0,ha,p0,pa,run_sarta);
  time1 = toc(tstart);

  tstart = tic;
  p2 = driver_sarta_cloud100layer_rtp(h0,ha,p0,pa,run_sarta);
  time2 = toc(tstart);

  run_sarta.ice_water_separator = -1;
  tstart = tic;
  p2X = driver_sarta_cloud100layer_rtp(h0,ha,p0,pa,run_sarta);
  time2X = toc(tstart);

  tstart = tic;
  pp = driver_pcrtm_cloud_rtp(h0,ha,p0,pa,run_sarta);
  timep = toc(tstart);

[time1 time2 time2X timep]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/CLOUD/
g = dogoodchan;
figure(1)
  plot(h0.vchan(g),rad2bt(h0.vchan(g),p0.robs1(g,:)) - rad2bt(h0.vchan(g),p0.rcalc(g,:)),'g',...
       h0.vchan(g),rad2bt(h0.vchan(g),p0.robs1(g,:)) - rad2bt(h0.vchan(g),p1.rcalc(g,:)),'b',...
       h0.vchan(g),rad2bt(h0.vchan(g),p0.robs1(g,:)) - rad2bt(h0.vchan(g),p2.rcalc(g,:)),'r',...
       h0.vchan(g),rad2bt(h0.vchan(g),p0.robs1(g,:)) - rad2bt(h0.vchan(g),p2X.rcalc(g,:)),'m',...
       h0.vchan(g),rad2bt(h0.vchan(g),p0.robs1(g,:)) - rad2bt(h0.vchan(g),pp.rcalc(g,:)),'k')

tobs = real(rad2bt(h0.vchan(g),p0.robs1(g,:)));    %% obs
t0   = rad2bt(h0.vchan(g),p0.rcalc(g,:));          %% orig cal with slabs
t1   = rad2bt(h0.vchan(g),p1.rcalc(g,:));          %% redoing the cal with slabs
t2   = rad2bt(h0.vchan(g),p2.rcalc(g,:));          %% srta100,klayers100 .. with ice_water_separator = +400
t2X   = rad2bt(h0.vchan(g),p2X.rcalc(g,:));        %% srta100,klayers100 .. with ice_water_separator = -1
tp   = rad2bt(h0.vchan(g),pp.rcalc(g,:));

figure(2)
  plot(h0.vchan(g),nanmean(tobs'-t0'),'g',h0.vchan(g),nanstd(tobs'-t0'),'g--',...
       h0.vchan(g),nanmean(tobs'-t1'),'b',h0.vchan(g),nanstd(tobs'-t1'),'b--',...
       h0.vchan(g),nanmean(tobs'-t2'),'r',h0.vchan(g),nanstd(tobs'-t2'),'r--',...
       h0.vchan(g),nanmean(tobs'-t2X'),'m',h0.vchan(g),nanstd(tobs'-t2X'),'m--',...
       h0.vchan(g),nanmean(tobs'-tp'),'k',h0.vchan(g),nanstd(tobs'-tp'),'k--')
axis([650 2700 -6 +6]); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tobs = real(rad2bt(h0.vchan,p0.robs1));
t0   = rad2bt(h0.vchan,p0.rcalc);
t1   = rad2bt(h0.vchan,p1.rcalc);
t2   = rad2bt(h0.vchan,p2.rcalc);
t2X  = rad2bt(h0.vchan,p2X.rcalc);
tp   = rad2bt(h0.vchan,pp.rcalc);

figure(3);
dbt = 200 : 1 : 300;
dbt = 200 : 2.5 : 300;
nobs = hist(tobs(903,:),dbt);
n0  = hist(t0(903,:),dbt);
n1  = hist(t1(903,:),dbt);
n2  = hist(t2(903,:),dbt);
n2X = hist(t2X(903,:),dbt);
np  = hist(tp(903,:),dbt);
  plot(dbt,nobs,'y','linewidth',4); hold on
  plot(dbt,n0,'g',dbt,n1,'b',dbt,n2,'r',n2X,dbt,'m',dbt,np,'k'); hold off
  hl = legend('Obs','Input','input rerun','sarta100 440mb','sarta100 no I/W separation',...
               'pcrtm100','location','northwest');
  set(hl,'fontsize',10); title('BT961')

figure(4);
dbt = 200 : 1 : 300;
dbt = 200 : 2.5 : 300;
dbt = 220 : 2.5 : 320;
nobs = hist(tobs(1291,:),dbt);
n0  = hist(t0(1291,:),dbt);
n1  = hist(t1(1291,:),dbt);
n2  = hist(t2(1291,:),dbt);
n2X = hist(t2X(1291,:),dbt);
np  = hist(tp(1291,:),dbt);
  plot(dbt,nobs,'y','linewidth',4); hold on
  plot(dbt,n0,'g',dbt,n1,'b',dbt,n2,'rx-',dbt,n2X,'mo-',dbt,np,'ks-'); hold off
  hl = legend('Obs','Input','input rerun','sarta100 440mb','sarta100 no I/W separation',...
               'pcrtm100','location','northwest');
  set(hl,'fontsize',10); title('BT1231')
