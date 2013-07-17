% testing the code
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil

iAIRSTestOneCld = +1;
iAIRSTestOneCld = -1;
 
if iAIRSTestOneCld > 0
  fname = '/asl/data/rtprod_airs/2012/05/01/cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp';
  [h,ha,p,pa] = rtpread(fname);
  [h,ha,p,pa] = rtpgrow(h,ha,p,pa);
else
  fname = '/strowdata1/shared/imbiriba/test_file_for_sergio.rtp';
  run_sarta.sartaclear_code = '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_wcon_nte';
  run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  [h,ha,p,pa] = rtpread(fname);
end

[h,p] = subset_rtp_allcloudfields(h,p,[],[],1:415);

run_sarta.cloud = +1;
run_sarta.clear = +1;
if iAIRSTestOneCld > 0
  [pall]  = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
  [px,ix] = driver_sarta_cloud_rtp_onecldtest(h,ha,p,pa,run_sarta,+1);  %% for water
  %[px,ix] = driver_sarta_cloud_rtp_onecldtest(h,ha,p,pa,run_sarta,-1);  %% for ice

  run_sarta.clear = +1;
  run_sarta.cloud = +1;
  run_sarta.ncol0 = -1;    %% so forces cfrac = 1
  addpath ../pcrtm
  addpath ../sarta
  cd ../pcrtm
  p1 = driver_pcrtm_cloud_rtp_onecldtest(h,ha,px,pa,run_sarta,+1);  %% for water
  %p1 = driver_pcrtm_cloud_rtp_onecldtest(h,ha,px,pa,run_sarta,-1);  %% for ice

  cd ../sarta
  p1.rsarta_water_cloud_only = px.rcalc;

  figure(1)
    plot(h.vchan,rad2bt(h.vchan,p1.rclearcalc)-rad2bt(h.vchan,p1.rcalc))
  figure(2)
    plot(h.vchan,rad2bt(h.vchan,p1.rclearcalc)-rad2bt(h.vchan,p1.rsarta_water_cloud_only))
    title('cloud effect')
  figure(3)
    %% rcalc is input, and NOT changed by driver_pcrtm_cloud_rtp
    plot(h.vchan,rad2bt(h.vchan,p1.rcalc)-rad2bt(h.vchan,p1.rsarta_water_cloud_only))  
    %% better be the same, both done by SARTA
    plot(h.vchan,rad2bt(h.vchan,p1.rclearcalc)-rad2bt(h.vchan,p1.sarta_clear))
    %% better be close, as this is clearsky
    plot(h.vchan,rad2bt(h.vchan,p1.rad_clrsky)-rad2bt(h.vchan,p1.sarta_clear))

    %% better be zero, as this is cloudy done by SARTA
    plot(h.vchan,rad2bt(h.vchan,p1.sarta_cloud)-rad2bt(h.vchan,p1.rsarta_water_cloud_only))

  figure(4)
    plot(h.vchan,rad2bt(h.vchan,p1.rad_allsky)-rad2bt(h.vchan,p1.sarta_cloud))
    title('PCRTM watercld - SARTA watercld')

  figure(2);
  figure(4);

  figure(3)
    boo = find(h.ichan == 1291);
    dbt = -10 : 1 : +30;
    n1 = hist(rad2bt(1231,p1.rclearcalc(boo,:))-rad2bt(1231,p1.sarta_cloud(boo,:)),dbt);
      n1 = n1/nansum(n1);
    n2 = hist(rad2bt(1231,p1.rad_allsky(boo,:))-rad2bt(1231,p1.sarta_cloud(boo,:)),dbt);
      n2 = n2/nansum(n2);
    n3 = hist(rad2bt(1231,pall.rcalc(boo,ix))-rad2bt(1231,px.rcalc(boo,:)),dbt);
      n3 = n3/nansum(n3);
    n4 = hist(rad2bt(1231,pall.rclearcalc(boo,ix))-rad2bt(1231,pall.rcalc(boo,ix)),dbt);
      n4 = n4/nansum(n4);
    plot(dbt,n1,'b',dbt,n2,'r',dbt,n3,'k',dbt,n4,'g')
    hl=legend('clear-watercld','PCRTM-SARTA','SARTAice/water-SARTAice','clear-ice/water');
    set(hl,'fontsize',8,'location','northwest');

    cloud_effect = sum(n1.*dbt);     cloud_effect_std = sqrt(sum(n1.*((dbt-cloud_effect).^2)));
    models_bias  = sum(n2.*dbt);     models_std       = sqrt(sum(n2.*((dbt-models_bias).^2)));
    fprintf(1,'cloud efect = %8.6f +/- %8.6f K \n',cloud_effect,cloud_effect_std);
    fprintf(1,'model diff  = %8.6f +/- %8.6f K \n',models_bias,models_std);

else
  [px] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
end
