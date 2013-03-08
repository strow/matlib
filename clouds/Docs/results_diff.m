chans6 = [359         445         903        1291        1557        1614];
chans38 = [        41          54         181         273         317         359         445         449 ...
                  532         758         903         904        1000        1020        1034        1055 ...
                 1075        1103        1249        1282        1291        1447        1475        1557 ...
                 1604        1614        1618        1660        1790        1866        1867        1868 ...
                 1878        1888        2112        2140        2321        2333];

clear p1 p2

p1.rad_allsky  = [];
p1.rad_clrsky  = [];
p1.sarta_cloud = [];
p1.sarta_clear = [];
p1.solzen = [];
p1.stemp = [];
p1.rlon = [];
p1.rlat = [];

p1.rad38_allsky  = [];
p1.rad38_clrsky  = [];
p1.sarta38_cloud = [];
p1.sarta38_clear = [];
p1.obs38         = [];
p1.solzen38 = [];
p1.stemp38 = [];
p1.rlon38 = [];
p1.rlat38 = [];

p1.pcrtm38_iceOD = [];
p1.pcrtm38_iceTOP = [];
p1.pcrtm38_iceDME = [];
p1.pcrtm38_waterOD = [];
p1.pcrtm38_waterTOP = [];
p1.pcrtm38_waterDME = [];

p1.sarta38_iceOD = [];
p1.sarta38_iceTOP = [];
p1.sarta38_iceDME = [];
p1.sarta38_waterOD = [];
p1.sarta38_waterTOP = [];
p1.sarta38_waterDME = [];

%%%%%%%%%%%%%%%%%%%%%%%%%

p2.rad_allsky  = [];
p2.rad_clrsky  = [];
p2.sarta_cloud = [];
p2.sarta_clear = [];
p2.solzen = [];
p2.stemp = [];
p2.rlon = [];
p2.rlat = [];

p2.rad38_allsky  = [];
p2.rad38_clrsky  = [];
p2.sarta38_cloud = [];
p2.sarta38_clear = [];
p2.obs38         = [];
p2.solzen38 = [];
p2.stemp38 = [];
p2.rlon38 = [];
p2.rlat38 = [];

p2.pcrtm38_iceOD = [];
p2.pcrtm38_iceTOP = [];
p2.pcrtm38_iceDME = [];
p2.pcrtm38_waterOD = [];
p2.pcrtm38_waterTOP = [];
p2.pcrtm38_waterDME = [];

p2.sarta38_iceOD = [];
p2.sarta38_iceTOP = [];
p2.sarta38_iceDME = [];
p2.sarta38_waterOD = [];
p2.sarta38_waterTOP = [];
p2.sarta38_waterDME = [];

%%%%%%%%%%%%%%%%%%%%%%%%%

ddStart = 1; ddEnd = 1;
for dd = ddStart : ddEnd
  fprintf(1,' day %2i \n',dd)
  fmain = ['/asl/data/rtprod_airs/2012/05/' num2str(dd,'%02d') '/cld_ecm_41ch.airs_ctr.2012.05.' num2str(dd,'%02d') '.'];
    for ii = 1 : 24
    foutA = [fmain  num2str(ii,'%02d') '_2378chansNEW.rtp'];
    foutB = [fmain  num2str(ii,'%02d') '_2378chansNEW_ncol0_1.rtp'];
    foutC = [fmain  num2str(ii,'%02d') '_2378chansNEW_ncol0_1_csum_0p3.rtp'];   %% fixed the random cld fracn when calling sarta

    fout1 = foutB;
    fout2 = foutC;

    ee1 = exist(fout1,'file');
    ee2 = exist(fout2,'file');
    if ee1 > 0 & ee2 > 0
      iYes = -1;

      [h,ha,p,pa] = oldrtpread(fout1);
      fprintf(1,'    %s \n',fout1);

      if length(p1.solzen) < 150000
        %% else sucks up too much memory
        iYes = 1;
        p1.rad_allsky  = [p1.rad_allsky   p.rad_allsky];
        p1.rad_clrsky  = [p1.rad_clrsky   p.rad_clrsky];
        p1.sarta_cloud = [p1.sarta_cloud  p.sarta_cloud];
        p1.sarta_clear = [p1.sarta_clear  p.sarta_clear];

        p1.solzen      = [p1.solzen       p.solzen];
        p1.stemp       = [p1.stemp        p.stemp];
        p1.rlon        = [p1.rlon         p.rlon];
        p1.rlat        = [p1.rlat         p.rlat];
      end

      p1.rad38_allsky  = [p1.rad38_allsky   p.rad_allsky(chans38,:)];
      p1.rad38_clrsky  = [p1.rad38_clrsky   p.rad_clrsky(chans38,:)];
      p1.sarta38_cloud = [p1.sarta38_cloud  p.sarta_cloud(chans38,:)];
      p1.sarta38_clear = [p1.sarta38_clear  p.sarta_clear(chans38,:)];
      p1.obs38         = [p1.obs38          p.robs1(chans38,:)];
      p1.solzen38      = [p1.solzen38       p.solzen];
      p1.stemp38       = [p1.stemp38        p.stemp];
      p1.rlon38        = [p1.rlon38         p.rlon];
      p1.rlat38        = [p1.rlat38         p.rlat];

      p1.pcrtm38_iceOD     = [p1.pcrtm38_iceOD         p.pcrtm_iceOD];
      p1.pcrtm38_iceTOP    = [p1.pcrtm38_iceTOP        p.pcrtm_iceCTOP];
      p1.pcrtm38_iceDME    = [p1.pcrtm38_iceDME        p.pcrtm_iceDME];
      p1.pcrtm38_waterOD   = [p1.pcrtm38_waterOD       p.pcrtm_waterOD];
      p1.pcrtm38_waterTOP  = [p1.pcrtm38_waterTOP      p.pcrtm_waterCTOP];
      p1.pcrtm38_waterDME  = [p1.pcrtm38_waterDME       p.pcrtm_waterDME];

      p1.sarta38_iceOD     = [p1.sarta38_iceOD         p.cngwat];
      p1.sarta38_iceTOP    = [p1.sarta38_iceTOP        p.cprtop];
      p1.sarta38_iceDME    = [p1.sarta38_iceDME        p.cpsize];
      p1.sarta38_waterOD   = [p1.sarta38_waterOD       p.cngwat2];
      p1.sarta38_waterTOP  = [p1.sarta38_waterTOP      p.cprtop2];
      p1.sarta38_waterDME  = [p1.sarta38_waterDME      p.cpsize2];

      [h,ha,p,pa] = oldrtpread(fout2);
      fprintf(1,'    %s \n\n',fout2);

      if iYes > 0
        %% else sucks up too much memory
        p2.rad_allsky  = [p2.rad_allsky   p.rad_allsky];
        p2.rad_clrsky  = [p2.rad_clrsky   p.rad_clrsky];
        p2.sarta_cloud = [p2.sarta_cloud  p.sarta_cloud];
        p2.sarta_clear = [p2.sarta_clear  p.sarta_clear];

        p2.solzen      = [p2.solzen       p.solzen];
        p2.stemp       = [p2.stemp        p.stemp];
        p2.rlon        = [p2.rlon         p.rlon];
        p2.rlat        = [p2.rlat         p.rlat];
      end

      p2.rad38_allsky  = [p2.rad38_allsky   p.rad_allsky(chans38,:)];
      p2.rad38_clrsky  = [p2.rad38_clrsky   p.rad_clrsky(chans38,:)];
      p2.sarta38_cloud = [p2.sarta38_cloud  p.sarta_cloud(chans38,:)];
      p2.sarta38_clear = [p2.sarta38_clear  p.sarta_clear(chans38,:)];
      p2.obs38         = [p2.obs38          p.robs1(chans38,:)];
      p2.solzen38      = [p2.solzen38       p.solzen];
      p2.stemp38       = [p2.stemp38        p.stemp];
      p2.rlon38        = [p2.rlon38         p.rlon];
      p2.rlat38        = [p2.rlat38         p.rlat];

      p2.pcrtm38_iceOD     = [p2.pcrtm38_iceOD         p.pcrtm_iceOD];
      p2.pcrtm38_iceTOP    = [p2.pcrtm38_iceTOP        p.pcrtm_iceCTOP];
      p2.pcrtm38_iceDME    = [p2.pcrtm38_iceDME        p.pcrtm_iceDME];
      p2.pcrtm38_waterOD   = [p2.pcrtm38_waterOD       p.pcrtm_waterOD];
      p2.pcrtm38_waterTOP  = [p2.pcrtm38_waterTOP      p.pcrtm_waterCTOP];
      p2.pcrtm38_waterDME  = [p2.pcrtm38_waterDME       p.pcrtm_waterDME];

      p2.sarta38_iceOD     = [p2.sarta38_iceOD         p.cngwat];
      p2.sarta38_iceTOP    = [p2.sarta38_iceTOP        p.cprtop];
      p2.sarta38_iceDME    = [p2.sarta38_iceDME        p.cpsize];
      p2.sarta38_waterOD   = [p2.sarta38_waterOD       p.cngwat2];
      p2.sarta38_waterTOP  = [p2.sarta38_waterTOP      p.cprtop2];
      p2.sarta38_waterDME  = [p2.sarta38_waterDME      p.cpsize2];

    end
  end
end

disp('done reading in data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /strowdata1/shared/sergio/MATLABCODE/SHOWSTATS/

fairs = instr_chans;
tS1 = rad2bt(fairs(chans38),p1.sarta38_cloud);  tS1 = real(tS1);
tP1 = rad2bt(fairs(chans38),p1.rad38_allsky);   tP1 = real(tP1);
tO1 = rad2bt(fairs(chans38),p1.obs38);          tO1 = real(tO1);

tS2 = rad2bt(fairs(chans38),p2.sarta38_cloud);  tS2 = real(tS2);
tP2 = rad2bt(fairs(chans38),p2.rad38_allsky);   tP2 = real(tP2);
if abs(nansum(nansum(p1.obs38-p2.obs38))) < 1e-3
  tO2 = tO1;
else
  tO2 = rad2bt(fairs(chans38),p2.obs38);          tO2 = real(tO2);
end

%% cind1 =         532         758         903        1249        1291        2321        2333
%% ff    =      0.8224      0.8999      0.9611      1.1295      1.2313      2.6037      2.6164

ii2616 = find(chans38 == 2333);
ii1231 = find(chans38 == 1291);
ii0961 = find(chans38 == 0903);
ii0820 = find(chans38 == 0532);

if length(p1.rlat) < 10000
  figure(1); plot(h.vchan,rad2bt(h.vchan,p1.sarta_cloud) - rad2bt(h.vchan,p1.rad_allsky)); title('SARTA - PCRTM allsky : v1')
  figure(2); plot(h.vchan,rad2bt(h.vchan,p1.sarta_clear) - rad2bt(h.vchan,p1.sarta_cloud)); title('SARTA : clr-cld : v1')

  figure(3); plot(h.vchan,rad2bt(h.vchan,p2.sarta_cloud) - rad2bt(h.vchan,p2.rad_allsky)); title('SARTA - PCRTM allskym : v2')
  figure(4); plot(h.vchan,rad2bt(h.vchan,p2.sarta_clear) - rad2bt(h.vchan,p2.sarta_cloud)); title('SARTA : clr-cld : v2')
else
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p1.sarta_cloud) - rad2bt(h.vchan,p1.rad_allsky),-20:1:+20,1,'SARTA - PCRTM allsky : v1');
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p1.sarta_clear) - rad2bt(h.vchan,p1.sarta_cloud),-10:1:+100,2,'SARTA : clr-cl : v1');

  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p2.sarta_cloud) - rad2bt(h.vchan,p2.rad_allsky),-20:1:+20,3,'SARTA - PCRTM allsky : v2');
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p2.sarta_clear) - rad2bt(h.vchan,p2.sarta_cloud),-10:1:+100,4,'SARTA : clr-cld : v2');
end

night = find(p1.solzen38 > 90);
day   = find(p1.solzen38 < 90);
oo = night;
addpath /asl/matlab2012/aslutil/

figure(5); simplemap(p1.rlat38(oo),p1.rlon38(oo),tO1(ii1231,oo));  title('v1 OBS BT 1231');
  caxis([200 320]); colorbar
figure(6); simplemap(p1.rlat38(oo),p1.rlon38(oo),tS1(ii1231,oo));  title('v1 SARTA cld BT 1231'); 
  caxis([200 320]); colorbar
figure(7); simplemap(p1.rlat38(oo),p1.rlon38(oo),tP1(ii1231,oo));  title('v1 PCRTM cld BT 1231'); 
  caxis([200 320]); colorbar
figure(8); simplemap(p2.rlat38(oo),p2.rlon38(oo),tS2(ii1231,oo));  title('v2 SARTA cld BT 1231'); 
  caxis([200 320]); colorbar
figure(9); simplemap(p2.rlat38(oo),p2.rlon38(oo),tP2(ii1231,oo));  title('v2 PCRTM cld BT 1231'); 
  caxis([200 320]); colorbar

dbt = -10 : 1 : +30;
tS = tS1; tP = tP1; 
  n2616dbt = hist(tS(ii2616,oo) - tP(ii2616,oo),dbt); n2616dbt = n2616dbt/sum(n2616dbt);
  n1231dbt = hist(tS(ii1231,oo) - tP(ii1231,oo),dbt); n1231dbt = n1231dbt/sum(n1231dbt);
  n0961dbt = hist(tS(ii0961,oo) - tP(ii0961,oo),dbt); n0961dbt = n0961dbt/sum(n0961dbt);
  n0820dbt = hist(tS(ii0820,oo) - tP(ii0820,oo),dbt); n0820dbt = n0820dbt/sum(n0820dbt);
figure(10); plot(dbt,[n0820dbt; n0961dbt; n1231dbt; n2616dbt],'linewidth',2); 
  title('v1 SARTA - PCRTM BTD'); hl = legend('820','960','1231','2616'); set(hl,'Fontsize',10)
  xlabel('dbt (K)'); ylabel('n(dbt)'); grid
tX = tS(ii2616,oo) - tP(ii2616,oo); 
tX = tS(ii1231,oo) - tP(ii1231,oo); 
nn = hist(tX,dbt); nn = nn/sum(nn); mu = trapz(dbt,dbt.*nn); sigma = sqrt(trapz(dbt,((dbt - mu).^2).*nn));
  [trapz(dbt,nn) mu sigma]

tS = tS2; tP = tP2; 
  n2616dbt = hist(tS(ii2616,oo) - tP(ii2616,oo),dbt); n2616dbt = n2616dbt/sum(n2616dbt);
  n1231dbt = hist(tS(ii1231,oo) - tP(ii1231,oo),dbt); n1231dbt = n1231dbt/sum(n1231dbt);
  n0961dbt = hist(tS(ii0961,oo) - tP(ii0961,oo),dbt); n0961dbt = n0961dbt/sum(n0961dbt);
  n0820dbt = hist(tS(ii0820,oo) - tP(ii0820,oo),dbt); n0820dbt = n0820dbt/sum(n0820dbt);
figure(11); plot(dbt,[n0820dbt; n0961dbt; n1231dbt; n2616dbt],'linewidth',2); 
  title('v2 SARTA - PCRTM BTD'); hl = legend('820','960','1231','2616'); set(hl,'Fontsize',10)
  xlabel('dbt (K)'); ylabel('n(dbt)'); grid
tX = tS(ii2616,oo) - tP(ii2616,oo); 
tX = tS(ii1231,oo) - tP(ii1231,oo); 
nn = hist(tX,dbt); nn = nn/sum(nn); mu = trapz(dbt,dbt.*nn); sigma = sqrt(trapz(dbt,((dbt - mu).^2).*nn));
  [trapz(dbt,nn) mu sigma]

tS = tS1; tP = tP1; 
figure(12); [ns rx ry] = myhist2d(tP(ii1231,oo),tS(ii1231,oo)-tP(ii1231,oo),200:1:320,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('PCRTM BT1231'); ylabel('dBT 1231 (K)'); text(323,32,'log10(cnt)'); title('v1 BT1231 : sarta - pcrtm')
tS = tS2; tP = tP2; 
figure(13); [ns rx ry] = myhist2d(tP(ii1231,oo),tS(ii1231,oo)-tP(ii1231,oo),200:1:320,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('PCRTM BT1231'); ylabel('dBT 1231 (K)'); text(323,32,'log10(cnt)'); title('v2 BT1231 : sarta - pcrtm')

tS = tS1; tP = tP1; 
figure(12); [ns rx ry] = myhist2d(tP(ii1231,oo),tS(ii1231,oo)-tP(ii1231,oo),200:1:320,-30:1:+30,+1); pcolor(rx,ry,(ns)); colorbar;
  xlabel('PCRTM BT1231 (K)'); ylabel('dBT 1231 (K)'); text(326,32,'cnt'); title('BT1231 : sarta - pcrtm')
  lala = colormap(jet); lala(1,:) = 1; colormap(lala); shading flat
  set(gca,'fontsize',10);
tS = tS2; tP = tP2; 
figure(13); [ns rx ry] = myhist2d(tP(ii1231,oo),tS(ii1231,oo)-tP(ii1231,oo),200:1:320,-30:1:+30,+1); pcolor(rx,ry,(ns)); colorbar;
  xlabel('PCRTM BT1231 (K)'); ylabel('dBT 1231 (K)'); text(326,32,'cnt'); title('BT1231 : sarta - pcrtm')
  lala = colormap(jet); lala(1,:) = 1; colormap(lala); shading flat
  set(gca,'fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%
tS = tS1; tP = tP1; xyz = p1.pcrtm38_iceOD; dxyz = 0:2:50;
figure(12); [ns rx ry] = myhist2d(xyz(oo),tS(ii1231,oo)-tP(ii1231,oo),dxyz,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('xyz'); ylabel('bias BT1231 : sarta-pcrtm'); text(323,32,'log10(cnt)'); title('v1 BT1231 : sarta - pcrtm')
tS = tS2; tP = tP2; xyz = p1.pcrtm38_iceOD; dxyz = 0:2:50;
figure(13); [ns rx ry] = myhist2d(xyz(oo),tS(ii1231,oo)-tP(ii1231,oo),dxyz,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('xyz'); ylabel('bias BT1231 : sarta-pcrtm'); text(323,32,'log10(cnt)'); title('v2 BT1231 : sarta - pcrtm')

tS = tS1; tP = tP1; xyz = p1.pcrtm38_iceOD - p1.sarta38_iceOD/50; dxyz = -25:1:25;
figure(12); [ns rx ry] = myhist2d(xyz(oo),tS(ii1231,oo)-tP(ii1231,oo),dxyz,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('xyz'); ylabel('bias BT1231 : sarta-pcrtm'); text(323,32,'log10(cnt)'); title('v1 BT1231 : sarta - pcrtm')
tS = tS2; tP = tP2; xyz = p1.pcrtm38_iceOD - p1.sarta38_iceOD/50; dxyz = -25:1:25;
figure(13); [ns rx ry] = myhist2d(xyz(oo),tS(ii1231,oo)-tP(ii1231,oo),dxyz,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('xyz'); ylabel('bias BT1231 : sarta-pcrtm'); text(323,32,'log10(cnt)'); title('v2 BT1231 : sarta - pcrtm')

tS = tS1; tP = tP1; xyz = p1.pcrtm38_iceTOP; dxyz = 100:50:500;
figure(12); [ns rx ry] = myhist2d(xyz(oo),tS(ii1231,oo)-tP(ii1231,oo),dxyz,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('xyz'); ylabel('bias BT1231 : sarta-pcrtm'); text(323,32,'log10(cnt)'); title('v1 BT1231 : sarta - pcrtm')
tS = tS2; tP = tP2; xyz = p2.pcrtm38_iceTOP; dxyz = 100:50:500;
figure(13); [ns rx ry] = myhist2d(xyz(oo),tS(ii1231,oo)-tP(ii1231,oo),dxyz,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('xyz'); ylabel('bias BT1231 : sarta-pcrtm'); text(323,32,'log10(cnt)'); title('v2 BT1231 : sarta - pcrtm')



site = transcom_match_sergio(p1.rlat38, p1.rlon38, 25);
oo = 1 : length(p1.rlat38);
oo = find(p1.solzen38 > 90);
oo = find((site == 13 | site == 14) & p1.solzen38 > 90);
tS = tS1; tP = tP1; tO = tO1;
figure(14);
  dbt = 180 : 1 : 330; nPCRTM = nanhist(tP(ii1231,oo),dbt); nPCRTM = nPCRTM/sum(nPCRTM);
  dbt = 180 : 1 : 330; nSARTA = nanhist(tS(ii1231,oo),dbt); nSARTA = nSARTA/sum(nSARTA);
  dbt = 180 : 1 : 330; nOBS   = nanhist(tO(ii1231,oo),dbt); nOBS   = nOBS/sum(nOBS);
  plot(dbt,nPCRTM,dbt,nSARTA,'r',dbt,nOBS,'k','linewidth',2); title('v1 : (b) PCRTM (r) SARTA (k) OBS'); grid
figure(15)
  plot(dbt,nOBS - nPCRTM,'b',dbt,nOBS - nSARTA,'r','linewidth',2); title('v1 : (b) OBS-PCRTM (r) OBS-SARTA')

site = transcom_match_sergio(p2.rlat38, p2.rlon38, 25);
oo = 1 : length(p2.rlat38);
oo = find(p2.solzen38 > 90);
oo = find((site == 13 | site == 14) & p2.solzen38 > 90);
tS = tS2; tP = tP2; tO = tO2;
figure(16);
  dbt = 180 : 1 : 330; nPCRTM = nanhist(tP(ii1231,oo),dbt); nPCRTM = nPCRTM/sum(nPCRTM);
  dbt = 180 : 1 : 330; nSARTA = nanhist(tS(ii1231,oo),dbt); nSARTA = nSARTA/sum(nSARTA);
  dbt = 180 : 1 : 330; nOBS   = nanhist(tO(ii1231,oo),dbt); nOBS   = nOBS/sum(nOBS);
  plot(dbt,nPCRTM,dbt,nSARTA,'r',dbt,nOBS,'k','linewidth',2); title('v2 : (b) PCRTM (r) SARTA (k) OBS'); grid
figure(17)
  plot(dbt,nOBS - nPCRTM,'b',dbt,nOBS - nSARTA,'r','linewidth',2); title('v2 : (b) OBS-PCRTM (r) OBS-SARTA')

iPrint = -1;
if iPrint > 0
  printfig(10,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm_calc_vs_sarta_calc_histV1','jpg')
  printfig(11,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm_calc_vs_sarta_calc_histV2','jpg')

  printfig(12,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm_calc_vs_sarta_biasV1','jpg')
  printfig(13,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm_calc_vs_sarta_biasV2','jpg')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
oo = find(p.ctype == 201); 
plot(p.pcrtm_iceDME(oo),p.cpsize(oo),'.',p.cpsize(oo),p.cpsize(oo),'k')
  xlabel('pcrtm ice deff');       ylabel('sarta ice deff');    
plot(p.pcrtm_iceCTOP(oo),p.cprtop(oo),'.',p.cprtop(oo),p.cprtop(oo),'k')
  xlabel('pcrtm ice ctop');       ylabel('sarta ice ctop');    
plot(p.pcrtm_iceOD(oo),p.cngwat(oo)/50,'.',p.cngwat(oo)/50,p.cngwat(oo)/50,'k')
  xlabel('pcrtm ice OD');       ylabel('sarta ice cngwat/50');    

oo = find(p.ctype == 101); 
plot(p.pcrtm_waterDME(oo),p.cpsize2(oo),'.',p.cpsize2(oo),p.cpsize2(oo),'k')
  xlabel('pcrtm water deff');       ylabel('sarta water deff');    
plot(p.pcrtm_waterCTOP(oo),p.cprtop2(oo),'.',p.cprtop2(oo),p.cprtop2(oo),'k')
  xlabel('pcrtm water ctop');       ylabel('sarta water ctop');    
plot(p.pcrtm_waterOD(oo),p.cngwat2(oo)/4,'.',p.cngwat2(oo)/4,p.cngwat2(oo)/4,'k')
  xlabel('pcrtm water OD');       ylabel('sarta water cngwat/4');    
%}
