%% see /strowdata1/s1/sergio/PCRTM_XIANGLEI/ECMWF2PCRTM_V2.1/sergio_make_nc.m

JOB = (datenum(2012,05,01,00,0,0));

[yy mm dd hh xjunk1 xjunk2] = datevec(JOB(1));

fprintf(1,'processing %4i %2i %2i : %2i \n',yy,mm,dd,hh)

addpath /strowdata1/shared/schou/prod_mat/gribtools  %(git repository)
addpath /asl/matlab/rtptools
addpath /asl/matlab/aslutil/
addpath /asl/matlab/science/
addpath /asl/matlab/h4tools/

thedateDIR = ['/asl/data/rtprod_airs/' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
dotstr = [num2str(yy) '.' num2str(mm,'%02d') '.' num2str(dd,'%02d') '.' num2str(hh,'%02d')];
thefilename = [thedateDIR '/cld_era_41ch.airs_ctr.' dotstr '.rtp'];

ncol0 = 50;
thefilenameOUT = [thedateDIR '/quick_test_cld_era_41ch.airs_ctr.' dotstr '.pcrtm.ncol' num2str(ncol0) '.rtp'];

ee = exist(thefilename);
eeOUT = exist(thefilenameOUT);
if eeOUT > 0
  rmer = ['!/bin/rm ' thefilenameOUT]; eval(rmer);
  eeOUT = exist(thefilenameOUT);
end

if ee == 0
  fprintf(1,'%10i\n',JOB)
  disp('rtp file DNE');
elseif eeOUT > 0
  fprintf(1,'%10i\n',JOB)
  disp('output file already exists');
else
  clear h ha p pa
  [h,ha,p0,pa] = rtpread(thefilename);
  p0 = simple_subset(50,p0);       %% for debug
  [h,ha,p,pa] = rtpgrow(h,ha,p0,pa);

  nboxes = length(p.stemp);  
  [nlev,nprof] = size(p.clwc);
  ncol = ncol0;
  overlap = 1;  %% maximum overlap
  overlap = 3;  %% maximum random overlap

  P   = double(p.plevs);  %% 1 mb = 1 hPa

  WCT = double(p.clwc);   %% cloud liquid water content in kg/kg
  ICT = double(p.ciwc);   %% cloud ice    water content in kg/kg
  cc  = double(p.cc);     %% cloud fraction

  TT = double(p.ptemp);        %% temperature profile
  if h.gunit(1) == 21 & h.gunit(2) == 21
    %% need to change from g/g to g/kg
    q  = 1000*double(p.gas_1);   %% wv profile in g/kg
    o3 = 1000*double(p.gas_3);   %% o3 profile in g/kg
  elseif h.gunit(1) == 20 & h.gunit(2) == 20
    %% already in g/kg
    q  = double(p.gas_1);   %% wv profile in g/kg
    o3 = double(p.gas_3);   %% o3 profile in g/kg
  else
    error('oops check h.gunit!!');
  end

  Ps = double(p.spres);   %% surface pressure in mb
  Ts = double(p.stemp);   %% surface temp in K
  
  sfctype = ones(size(p.stemp)) * -9999;   %% so that we can specify emissivity
  efreq   = double(p.efreq);
  emis    = double(p.emis);

  zen_ang = double(p.scanang);
  co2     = ones(size(p.stemp)) * (370 + (yy-2002)*2.2);

  % use_Xiuhong  %% for debug default 2012/05/01  00:00-01:00 UTC
 
  yymmddhhstr = datestr(JOB(1),'yyyymmddHH');
  parname = ['/strowdata1/s1/sergio/PCRTM_XIANGLEI/NEWVERS/PCRTM2AIRS_spec/JUNK/pcrtm' yymmddhhstr '.tmp'];
  ppath = '/strowdata1/s1/sergio/PCRTM_XIANGLEI/NEWVERS/PCRTM_V2.1_for_AIRS/code_changed/Run/';
  parnameout = [parname '.out'];

  whos P WCT ICE cc TT q o3 Ps Ts sfctype efreq emis zen_ang co2 
  fprintf(1,'making PCRTM input file %s \n',parname)

  [rad_allsky rad_clrsky tmpjunk] = PCRTM_compute_for_AIRS_spectra(nboxes,nlev, ncol, overlap, ...
                                                           P, WCT, ICT, cc, TT, q, o3, Ps, Ts, ...
                                                           sfctype,efreq,emis, ...
                                                           zen_ang,co2, parname,ppath);

  get_sarta_clear
  get_extra_p0_params
  rtpwrite(thefilenameOUT,h,ha,p0,pa);

%  rmer = ['!/bin/rm ' parname ' ' parnameout];
%  eval(rmer);

  figure(1); plot(h.vchan,rad2bt(h.vchan,p0.rcalc)); title('SARTA CLD');
  figure(2); plot(h.vchan,rad2bt(h.vchan,p0.rad_allsky)); title('PCRTM allsky');
  figure(3); plot(h.vchan,rad2bt(h.vchan,p0.rad_clrsky)); title('PCRTM clrsky');

  figure(4); plot(h.vchan,rad2bt(h.vchan,p0.rad_clrsky) - rad2bt(h.vchan,p0.rad_allsky)); 
    title('PCRTM : clr - cld')
  figure(5); plot(h.vchan,rad2bt(h.vchan,p0.sarta_clear) - rad2bt(h.vchan,p0.rcalc)); 
    title('SARTA : clr - cld')

  figure(6); plot(h.vchan,rad2bt(h.vchan,p0.rad_clrsky) - rad2bt(h.vchan,p0.sarta_clear)); 
    title('clr : PCRTM - SARTA')
  figure(7); plot(h.vchan,rad2bt(h.vchan,p0.rad_allsky) - rad2bt(h.vchan,p0.rcalc)); 
    title('cld : PCRTM  - SARTA')

  aa = load('xianglei_2012_05_01_00hrs.mat');

end
disp('LOOK AT JUNK/pcrtm2012050100.tmp')

display_quick_test_driver
