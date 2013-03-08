chans6 = [359         445         903        1291        1557        1614];
chans38 = [        41          54         181         273         317         359         445         449 ...
                  532         758         903         904        1000        1020        1034        1055 ...
                 1075        1103        1249        1282        1291        1447        1475        1557 ...
                 1604        1614        1618        1660        1790        1866        1867        1868 ...
                 1878        1888        2112        2140        2321        2333];

clear p1

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


for dd = 1 : 31
  fprintf(1,' day %2i \n',dd)
  fmain = ['/asl/data/rtprod_airs/2012/05/' num2str(dd,'%02d') '/cld_ecm_41ch.airs_ctr.2012.05.' num2str(dd,'%02d') '.'];
    for ii = 1 : 24
    fout = [fmain  num2str(ii,'%02d') '_2378chans.rtp'];
    fout = [fmain  num2str(ii,'%02d') '_2378chansNEW.rtp'];
    ee = exist(fout,'file');
    if ee > 0
      [h,ha,p,pa] = oldrtpread(fout);
      fprintf(1,'    %s \n',fout);

      if length(p1.solzen) < 150000
        %% else sucks up too much memory
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

    end
  end
end

addpath /strowdata1/shared/sergio/MATLABCODE/SHOWSTATS/

if length(p1.rlat) < 10000
  figure(1); plot(h.vchan,rad2bt(h.vchan,p1.sarta_clear)); title('SARTA clear')
  figure(2); plot(h.vchan,rad2bt(h.vchan,p1.sarta_clear) - rad2bt(h.vchan,p1.rad_clrsky)); title('SARTA - PCRTM clear')
  figure(3); plot(h.vchan,rad2bt(h.vchan,p1.sarta_cloud)); title('SARTA cloud')
  figure(4); plot(h.vchan,rad2bt(h.vchan,p1.sarta_cloud) - rad2bt(h.vchan,p1.rad_allsky)); title('SARTA - PCRTM allsky')
  figure(5); plot(h.vchan,rad2bt(h.vchan,p1.sarta_clear) - rad2bt(h.vchan,p1.sarta_cloud)); title('SARTA : clr-cld')
else
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p1.sarta_clear),200:1:320,1,'SARTA clear'); 
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p1.sarta_clear) - rad2bt(h.vchan,p1.rad_clrsky),-2.5:0.25:+2.5,2,'SARTA - PCRTM clear');
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p1.sarta_cloud),200:1:320,3,'SARTA cloud');
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p1.sarta_cloud) - rad2bt(h.vchan,p1.rad_allsky),-20:1:+20,4,'SARTA - PCRTM allsky');
  na = show2dstatsVers2(h.vchan,rad2bt(h.vchan,p1.sarta_clear) - rad2bt(h.vchan,p1.sarta_cloud),-10:1:+100,5,'SARTA : clr-cld');
end

fairs = instr_chans;
tS = rad2bt(fairs(chans38),p1.sarta38_cloud);  tS = real(tS);
tP = rad2bt(fairs(chans38),p1.rad38_allsky);   tP = real(tP);
tO = rad2bt(fairs(chans38),p1.obs38);          tO = real(tO);


%% cind1 =         532         758         903        1249        1291        2321        2333
%% ff    =      0.8224      0.8999      0.9611      1.1295      1.2313      2.6037      2.6164

ii2616 = find(chans38 == 2333);
ii1231 = find(chans38 == 1291);
ii0961 = find(chans38 == 0903);
ii0820 = find(chans38 == 0532);

night = find(p1.solzen38 > 90);
day   = find(p1.solzen38 < 90);
oo = night;
addpath /asl/matlab2012/aslutil/
figure(6); simplemap(p1.rlat38(oo),p1.rlon38(oo),tO(ii1231,oo));  title('OBS BT 1231');
  caxis([200 320]); colorbar
figure(7); simplemap(p1.rlat38(oo),p1.rlon38(oo),tS(ii1231,oo));  title('SARTA cld BT 1231'); 
  caxis([200 320]); colorbar
figure(8); simplemap(p1.rlat38(oo),p1.rlon38(oo),tP(ii1231,oo));  title('PCRTM cld BT 1231'); 
  caxis([200 320]); colorbar
figure(9); simplemap(p1.rlat38(oo),p1.rlon38(oo),tS(ii1231,oo)-tP(ii1231,oo));  title('SARTA - PCRTM cld BT 1231'); 
  caxis([-5 +5]); colorbar;

dbt = -10 : 1 : +30; 
  n2616dbt = hist(tS(ii2616,oo) - tP(ii2616,oo),dbt); n2616dbt = n2616dbt/sum(n2616dbt);
  n1231dbt = hist(tS(ii1231,oo) - tP(ii1231,oo),dbt); n1231dbt = n1231dbt/sum(n1231dbt);
  n0961dbt = hist(tS(ii0961,oo) - tP(ii0961,oo),dbt); n0961dbt = n0961dbt/sum(n0961dbt);
  n0820dbt = hist(tS(ii0820,oo) - tP(ii0820,oo),dbt); n0820dbt = n0820dbt/sum(n0820dbt);
figure(10); plot(dbt,[n0820dbt; n0961dbt; n1231dbt; n2616dbt],'linewidth',2); 
  title('SARTA - PCRTM BTD'); hl = legend('820','960','1231','2616'); set(hl,'Fontsize',10)
  xlabel('dbt (K)'); ylabel('n(dbt)'); grid
tX = tS(ii2616,oo) - tP(ii2616,oo); 
tX = tS(ii1231,oo) - tP(ii1231,oo); 

nn = hist(tX,dbt); nn = nn/sum(nn); mu = trapz(dbt,dbt.*nn); sigma = sqrt(trapz(dbt,((dbt - mu).^2).*nn));
  [trapz(dbt,nn) mu sigma]


figure(11); plot(tP(ii1231,oo),tP(ii1231,oo)-tS(ii1231,oo),'.');
figure(12); [ns rx ry] = myhist2d(tP(ii1231,oo),tS(ii1231,oo)-tP(ii1231,oo),200:1:320,-30:1:+30,-1); pcolor(rx,ry,log10(ns)); colorbar;
  xlabel('PCRTM BT1231'); ylabel('bias BT1231 : sarta-pcrtm'); text(323,32,'log10(cnt)'); title('BT1231 : sarta - pcrtm')

site = transcom_match_sergio(p1.rlat38, p1.rlon38, 25);
oo = 1 : length(p1.rlat38);
oo = find(p1.solzen38 > 90);
oo = find((site == 13 | site == 14) & p1.solzen38 > 90);
figure(13);
  dbt = 180 : 1 : 330; nPCRTM = nanhist(tP(ii1231,oo),dbt); nPCRTM = nPCRTM/sum(nPCRTM);
  dbt = 180 : 1 : 330; nSARTA = nanhist(tS(ii1231,oo),dbt); nSARTA = nSARTA/sum(nSARTA);
  dbt = 180 : 1 : 330; nOBS   = nanhist(tO(ii1231,oo),dbt); nOBS   = nOBS/sum(nOBS);
  plot(dbt,nPCRTM,dbt,nSARTA,'r',dbt,nOBS,'k','linewidth',2); title('(b) PCRTM (r) SARTA (k) OBS'); grid
figure(14)
  plot(dbt,nOBS - nPCRTM,'b',dbt,nOBS - nSARTA,'r','linewidth',2); title('(b) OBS-PCRTM (r) OBS-SARTA')

iPrint = -1;
if iPrint > 0
  printfig(4,'~/MATLABCODE/matlib/clouds/Docs/Figs/spectra_sartaVSpcrtm','jpg')
  printfig(6,'~/MATLABCODE/matlib/clouds/Docs/Figs/obs1231_simplemap','jpg')
  printfig(7,'~/MATLABCODE/matlib/clouds/Docs/Figs/sarta1231_simplemap','jpg')
  printfig(8,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm1231_simplemap','jpg')
  printfig(10,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm_calc_vs_sarta_calc_hist','jpg')
  printfig(12,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm_calc_vs_sarta_bias','jpg')
  printfig(13,'~/MATLABCODE/matlib/clouds/Docs/Figs/pcrtm_vs_sarta_vs_obs_pdf1231','jpg')
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