robs      = [];
yrobs      = [];

rclrSARTA = [];
rclrPCRTM = [];
rcldSARTA = [];
rcldPCRTM = [];

yrclrSARTA = [];
yrclrPCRTM = [];
yrcldSARTA = [];
yrcldPCRTM = [];

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

clear t* oo gg

dd = 01; hh = 00;
iERAorECMWF = input('Enter (+1) ERA (-1) ECM : ');

fname0 = ['/asl/data/rtprod_airs/2012/05/' num2str(dd,'%02d')];
if iERAorECMWF == -1
  fname = [fname0 '/cld_ecm_41ch.airs_ctr.2012.05.'];
  fname = [fname0 '/xcld_ecm_41ch.airs_ctr.2012.05.'];
  fnamey = [fname0 '/ycld_ecm_41ch.airs_ctr.2012.05.'];
elseif iERAorECMWF == +1
  fname = [fname0 '/cld_era_41ch.airs_ctr.2012.05.'];
  fname = [fname0 '/xcld_era_41ch.airs_ctr.2012.05.'];
  fnamey = [fname0 '/ycld_era_41ch.airs_ctr.2012.05.'];
end
fname  = [fname  num2str(dd,'%02d') '.' num2str(hh,'%02d') '.pcrtm.ncol50.rtp'];
fnamey = [fnamey  num2str(dd,'%02d') '.' num2str(hh,'%02d') '.pcrtm.ncol50.rtp'];
  [h,ha,p,pa]   = rtpread(fname);
  [hy,ha,py,pa] = rtpread(fnamey);

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

  bonk = [];
  ybonk = [];

  for hh = 0 : 24
    fname0 = ['/asl/data/rtprod_airs/2012/05/' num2str(dd,'%02d')];
    if iERAorECMWF == -1
      fname = [fname0 '/cld_ecm_41ch.airs_ctr.2012.05.'];
      fname = [fname0 '/xcld_ecm_41ch.airs_ctr.2012.05.'];
      fnamey = [fname0 '/xcld_ecm_41ch.airs_ctr.2012.05.'];
    elseif iERAorECMWF == +1
      fname = [fname0 '/cld_era_41ch.airs_ctr.2012.05.'];
      fname = [fname0 '/xcld_era_41ch.airs_ctr.2012.05.'];
      fnamey = [fname0 '/ycld_era_41ch.airs_ctr.2012.05.'];
    end
    fname  = [fname  num2str(dd,'%02d') '.' num2str(hh,'%02d') '.pcrtm.ncol50.rtp'];
    fnamey = [fnamey  num2str(dd,'%02d') '.' num2str(hh,'%02d') '.pcrtm.ncol50.rtp'];
    ee = exist(fname);
    eey = exist(fnamey);
    if ee > 0 & eey > 0
      [h,ha,p,pa] = rtpread(fname);
      [hy,ha,py,pa] = rtpread(fnamey);

      find_ix

      if length(ix) > 0
        robs      = [robs      p.robs1(oo,ix)];
        rclrSARTA = [rclrSARTA p.sarta_clear(oo,ix)];
        rcldSARTA = [rcldSARTA p.rcalc(oo,ix)];

        rclrPCRTM = [rclrPCRTM p.rad_clrsky(oo,ix)];
        rcldPCRTM = [rcldPCRTM p.rad_allsky(oo,ix)];
        bonk      = [bonk p.rad_allsky(oo,ix)];

        yrclrPCRTM = [yrclrPCRTM py.rad_clrsky(oo,ix)];
        yrcldPCRTM = [yrcldPCRTM py.rad_allsky(oo,ix)];
        ybonk      = [ybonk py.rad_allsky(oo,ix)];

        nonk1 = nansum(nansum(p.robs1(oo,ix) - py.robs1(oo,ix)));
        nonk2 = nansum(nansum(p.rcalc(oo,ix) - py.rcalc(oo,ix)));
        nonk3 = nansum(nansum(p.sarta_clear(oo,ix) - py.sarta_clear(oo,ix)));
        if (abs(nonk1) > eps) | (abs(nonk2) > eps) | (abs(nonk3) > eps)
          fprintf(1,'%s fname has different x,y data \n',fname);
        end

boo = find(isnan(p.rad_allsky(oo,ix)));
if length(boo) > 0
  fprintf(1,'oops %s has %4i NaNs for PCRTM1231 allsky \n',fname,length(boo))
  plot(p.rad_allsky(oo,ix))
end

        lat       = [lat p.rlat(ix)];
        lon       = [lon p.rlon(ix)];
        landfrac  = [landfrac p.landfrac(ix)];
        stemp     = [stemp p.stemp(ix)];

        cloud_stuff

      end   %% if length(ix) > 0
    end     %% if ee 
  end       %% for hh
  datx = real(rad2bt(h.vchan(oo),bonk));
  daty = real(rad2bt(h.vchan(oo),ybonk));
  figure(1); plot(datx-daty); 
  figure(2); plot(datx); 
  disp('ret to continue'); pause
end         %% for dd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE/SHOWSTATS

fname 
fnamey

datx = real(rad2bt(h.vchan(oo),rcldPCRTM));
daty = real(rad2bt(h.vchan(oo),yrcldPCRTM));

[n x y]=myhist2d(datx,daty,[200:1:320],[200:1:320],-1);    
pcolor(x,y,n); shading flat; title('BT 1231');
gaga = colormap; gaga(1,:) = 1; colormap(gaga); caxis([0 50]); colorbar('horiz');
xlabel('orig PCRTM','Fontsize',12); ylabel('new PCRTM','Fontsize',12); 

