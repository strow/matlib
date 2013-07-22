% testing the code
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath ../sarta

run_sarta.clear = +1;
run_sarta.cloud = +1;
[h,ha,p,pa] = rtpread('/asl/data/rtprod_airs/2012/05/01/cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp');
[h,ha,p,pa] = rtpgrow(h,ha,p,pa);
[h,p] = subset_rtp_allcloudfields(h,p,[],[],10);

run_sarta.ncol0 = -1;
tic
p1 = driver_pcrtm_cloud_rtp(h,ha,p,pa,run_sarta);
toc
rtpwrite('/asl/data/rtprod_airs/2012/05/01/pcrtm_cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp',h,ha,p1,pa);

i1231 = abs(h.vchan - 1231); i1231 = find(i1231 == min(i1231));
figure(1); 
  scatter_coast(p.rlon,p.rlat,10,rad2bt(h.vchan(i1231),p1.sarta_clear(i1231,:))); 
  title('SARTA clear for 1231 cm-1')
figure(2); 
  scatter_coast(p.rlon,p.rlat,10,rad2bt(h.vchan(i1231),p1.rad_clrsky(i1231,:))); 
  title('PCRTM clear for 1231 cm-1')
figure(3); 
  scatter_coast(p.rlon,p.rlat,10,rad2bt(h.vchan(i1231),p1.rad_allsky(i1231,:))); 
  title('PCRTM cloud for 1231 cm-1')

figure(1); 
  scatter_coast(p.rlon,p.rlat,10,rad2bt(h.vchan(i1231),p1.sarta_clear(i1231,:))); 
  title('SARTA clear for 1231 cm-1')
figure(2); 
  scatter_coast(p.rlon,p.rlat,10,...
      rad2bt(h.vchan(i1231),p1.sarta_clear(i1231,:))-rad2bt(h.vchan(i1231),p1.rad_clrsky(i1231,:)))
  title('clear : SARTA-PCRTM for 1231 cm-1')
figure(3); 
  scatter_coast(p.rlon,p.rlat,10,...
      rad2bt(h.vchan(i1231),p1.sarta_clear(i1231,:))-rad2bt(h.vchan(i1231),p1.rad_allsky(i1231,:)))
  title('cld forcing : SARTAclr-PCRTMcld for 1231 cm-1')

figure(1)
  plot(h.vchan,rad2bt(h.vchan,p1.sarta_clear),'b',h.vchan,rad2bt(h.vchan,p1.rad_clrsky),'r')
  title('CLR SKY : (b) kCARTA (r) PCRTM')
  plot(h.vchan,rad2bt(h.vchan,p1.sarta_clear)-rad2bt(h.vchan,p1.rad_clrsky))
  title('CLR SKY : kCARTA - PCRTM')
  dBT = -2 : 0.1 : +2;
  clear nn
  for ii = 1 : h.nchan
    z = rad2bt(h.vchan(ii),p1.sarta_clear(ii,:))-rad2bt(h.vchan(ii),p1.rad_clrsky(ii,:));
    nn(ii,:) = hist(z,dBT);
  end
  pcolor(h.vchan,dBT,nn'); shading flat; xlabel('freq cm-1'); ylabel('bias K'); colorbar
