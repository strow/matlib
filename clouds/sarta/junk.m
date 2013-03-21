  addpath /asl/matlib/h4tools
  addpath /asl/matlib/rtptools
  addpath /asl/matlib/aslutil

  iAIRS = +1;
  iAIRS = -1;
 
  if iAIRS > 0
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
  if iAIRS > 0
    [px,ix] = driver_sarta_cloud_rtp_watercldtest(h,ha,p,pa,run_sarta);

    run_sarta.clear = +1;
    run_sarta.cloud = +1;
    addpath ../pcrtm
    addpath ../sarta
    cd ../pcrtm
    p1 = driver_pcrtm_cloud_rtp(h,ha,px,pa,run_sarta);
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

  else
    [px] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
  end
