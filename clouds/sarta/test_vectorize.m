%[pse0,orig_slabs] = driver_sarta_cloud_rtp(h,hx,profsub,pattr,run_sarta);   %% code before vectorization
%[pse,orig_slabs]  = driver_sarta_cloud_rtp(h,hx,profsub,pattr,run_sarta);   %% code after vectorization
%[psex0,orig_slabs] = driver_sarta_cloud_rtp(h,hx,profsub,pattr,run_sarta);  %% trying to go back to code before vectorization

figure(1)
clf
subplot(221)
  plot(pse0.cpsize,pse.cpsize,'b.',pse0.cpsize2,pse.cpsize2,'r.',pse0.cpsize,pse0x.cpsize,'c.',pse0.cpsize2,pse0x.cpsize2,'m.'); title('cpsize')
subplot(222)
  plot(pse0.cfrac,pse.cfrac,'b.',pse0.cfrac2,pse.cfrac2,'r.',pse0.cfrac,pse0x.cfrac,'c.',pse0.cfrac2,pse0x.cfrac2,'m.'); title('cfrac')
subplot(223)
  plot(pse0.cngwat,pse.cngwat,'b.',pse0.cngwat2,pse.cngwat2,'r.',pse0.cngwat,pse0x.cngwat,'c.',pse0.cngwat2,pse0x.cngwat2,'m.'); title('cngwat')
subplot(224)
  plot(pse0.cprtop,pse.cprtop,'b.',pse0.cprtop2,pse.cprtop2,'r.',pse0.cprtop,pse0x.cprtop,'c.',pse0.cprtop2,pse0x.cprtop2,'m.'); title('cprtop')
  axis([0 1000 0 1000])

figure(2); clf
clf
subplot(221)
  plot(pse0.cpsize,psex0.cpsize,'b.',pse0.cpsize2,psex0.cpsize2,'r.',pse0.cpsize,pse0x.cpsize,'c.',pse0.cpsize2,pse0x.cpsize2,'m.'); title('cpsize')
subplot(222)
  plot(pse0.cfrac,psex0.cfrac,'b.',pse0.cfrac2,psex0.cfrac2,'r.',pse0.cfrac,pse0x.cfrac,'c.',pse0.cfrac2,pse0x.cfrac2,'m.'); title('cfrac')
subplot(223)
  plot(pse0.cngwat,psex0.cngwat,'b.',pse0.cngwat2,psex0.cngwat2,'r.',pse0.cngwat,pse0x.cngwat,'c.',pse0.cngwat2,pse0x.cngwat2,'m.'); title('cngwat')
subplot(224)
  plot(pse0.cprtop,psex0.cprtop,'b.',pse0.cprtop2,psex0.cprtop2,'r.',pse0.cprtop,pse0x.cprtop,'c.',pse0.cprtop2,pse0x.cprtop2,'m.'); title('cprtop')
  axis([0 1000 0 1000])

figure(3); clf
plot(pse0.cngwat,pse0.cngwat-psex0.cngwat,'.'); grid
ijk = find(pse0.cngwat > 5 & pse0.cngwat-psex0.cngwat > 2)
[hwierd,pwierd] = subset_rtp_allcloudfields(h,profsub,[],[],ijk);

%% keep changing iWhichInterp
run_sarta.iWhichInterp = 0;  %% intepr1
  [pw0,orig_slabsw] = driver_sarta_cloud_rtp(hwierd,hx,pwierd,pattr,run_sarta);
run_sarta.iWhichInterp = 1;  %% interp1qr
  [pw1,orig_slabsw] = driver_sarta_cloud_rtp(hwierd,hx,pwierd,pattr,run_sarta);

disp('ctype    cngwat    cprtop    ctype2     cngwat2    cprtop2');
disp('----------------------------------------------------------')
[pw0.ctype    pw1.ctype    psex0.ctype(ijk)    pse0.ctype(ijk); ...
 pw0.cngwat   pw1.cngwat   psex0.cngwat(ijk)   pse0.cngwat(ijk); ...
 pw0.cprtop   pw1.cprtop   psex0.cprtop(ijk)   pse0.cprtop(ijk); ...
 pw0.ctype2   pw1.ctype2   psex0.ctype2(ijk)   pse0.ctype2(ijk);...
 pw0.cngwat2  pw1.cngwat2  psex0.cngwat2(ijk)  pse0.cngwat2(ijk);...
 pw0.cprtop2  pw1.cprtop2  psex0.cprtop2(ijk)  pse0.cprtop2(ijk); ...
]'

figure(4); clf
plot(pse0.cngwat,pse0.cngwat-psex0.cngwat,'.'); grid
ijk = find(abs(pse0.cngwat-psex0.cngwat) > 0.1); whos ijk
[hwierd,pwierd] = subset_rtp_allcloudfields(h,profsub,[],[],ijk);
%% keep changing iWhichInterp
run_sarta.iWhichInterp = 0;  %% intepr1
  [pw0,orig_slabsw] = driver_sarta_cloud_rtp(hwierd,hx,pwierd,pattr,run_sarta);
run_sarta.iWhichInterp = 1;  %% interp1qr
  [pw1,orig_slabsw] = driver_sarta_cloud_rtp(hwierd,hx,pwierd,pattr,run_sarta);
plot(pw0.cngwat,pw0.cngwat-pw1.cngwat,'.',pw0.cngwat,pw0.cngwat-pse0.cngwat(ijk),'r.')
plot(pw0.cngwat,pw0.cngwat-pw1.cngwat2,'.',pw0.cngwat,pw0.cngwat-pw1.cngwat,'r.')
