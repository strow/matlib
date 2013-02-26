%% this driver file can run off calculations using Xianglei's codes

%% need the "x" in there because Paul's clustcmd assumes something for numbers of the form YYYYMMDDHH
%% so need to fool it a little
%%
% local running to test
% clustcmd -L xdriver_PCRTM_compute_for_AIRS_spectra_ERAorECMWF.m YYYYMMDDHH0:YYYYMMDDHHF
% eg 2012050100
% clustcmd -L xdriver_PCRTM_compute_for_AIRS_spectra_ERAorECMWF.m YYYMMDDHH0:YYYYMMDDHHF
%
% otherwise when happy
% clustcmd -q medium -n 12 -p 1 xdriver_PCRTM_compute_for_AIRS_spectra_ERAorECMWF.m 2012050100:2012053123
%
% or
% clustcmd -q long_contrib -n 12 -p 1 xdriver_PCRTM_compute_for_AIRS_spectra_ERAorECMWF.m 2012050100:2012053123
% can test with JOB = datenum(2012,06,01,00,0,0);
% [yy mm dd hh] = datevec(datenum(2012,06,01,00,0,0));
%

%yymmddhhstr = num2str(JOB);
%yy = str2num(yymmddhhstr(01:04));
%mm = str2num(yymmddhhstr(05:06));
%dd = str2num(yymmddhhstr(07:08));
%hh = str2num(yymmddhhstr(09:10));
%if hh > 23
%  disp('can only do jobs with hh between 00 and 23');
%  return
%end

[yy mm dd hh xjunk1 xjunk2] = datevec(JOB(1));

fprintf(1,'processing %4i %2i %2i : %2i \n',yy,mm,dd,hh)

addpath /strowdata1/shared/schou/prod_mat/gribtools  %(git repository)
addpath /asl/matlab/rtptools
addpath /asl/matlab/aslutil/
addpath /asl/matlab/science/
addpath /asl/matlab/h4tools/

iERAorECMWF = -1; %% do ECM
iERAorECMWF = +1; %% do ERA

thedateDIR = ['/asl/data/rtprod_airs/' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
dotstr = [num2str(yy) '.' num2str(mm,'%02d') '.' num2str(dd,'%02d') '.' num2str(hh,'%02d')];
if iERAorECMWF == 1
  thefilename = [thedateDIR '/cld_era_41ch.airs_ctr.' dotstr '.rtp'];
elseif iERAorECMWF == -1
  thefilename = [thedateDIR '/cld_ecm_41ch.airs_ctr.' dotstr '.rtp'];
end

ncol0 = 50;
if iERAorECMWF == 1
  thefilenameOUT = [thedateDIR '/xcld_era_41ch.airs_ctr.']; %% mine
  thefilenameOUT = [thedateDIR '/ycld_era_41ch.airs_ctr.']; %% xiuhong
elseif iERAorECMWF == -1
  thefilenameOUT = [thedateDIR '/xcld_ecm_41ch.airs_ctr.']; %% mine
  thefilenameOUT = [thedateDIR '/ycld_ecm_41ch.airs_ctr.']; %% xiuhong
end

thefilenameOUT = [thefilenameOUT  dotstr '.pcrtm.ncol' num2str(ncol0) '.rtp'];

iChunk = 100;

ee = exist(thefilename);
eeOUT = exist(thefilenameOUT);
if ee == 0
  fprintf(1,'%10i\n',JOB)
  disp('rtp file DNE');
elseif eeOUT > 0
  fprintf(1,'%10i %s\n',JOB,thefilenameOUT)
  disp('output file already exists');
else
  clear h ha p pa
  [h,ha,p0ALL,pa] = rtpread(thefilename);
  p0ALLX = p0ALL;

  for iInd = 1 : ceil(length(p0ALL.xtrack)/iChunk);

    inds = (1:iChunk) + (iInd-1)*iChunk;
    inds = intersect(1:length(p0ALL.xtrack),inds);
    p0 = index_subset(inds,p0ALLX); 

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
    fprintf(1,'making PCRTM input file %s for iChunk %3i \n',parname,iInd)

    [rad_allsky rad_clrsky tmpjunk] = PCRTM_compute_for_AIRS_spectra(nboxes,nlev, ncol, overlap, ...
                                                           P, WCT, ICT, cc, TT, q, o3, Ps, Ts, ...
                                                           sfctype,efreq,emis, ...
                                                           zen_ang,co2, parname,ppath);
%    [rad_allsky2 rad_clrsky tmpjunk2] = PCRTM_compute_for_AIRS_spectra_V2(nboxes,nlev, ncol, overlap, ...
%                                                           P, WCT, ICT, cc, TT, q, o3, Ps, Ts, ...
%                                                           sfctype,efreq,emis, ...
%                                                           zen_ang,co2, parname,ppath);

    get_sarta_clear
    get_extra_p0_params_inds
    rmer = ['!/bin/rm ' parname ' ' parnameout];
    eval(rmer);

  end

  rtpwrite(thefilenameOUT,h,ha,p0ALL,pa);

  figure(1); plot(h.vchan,rad2bt(h.vchan,p0ALL.rcalc)); title('SARTA CLD');
  figure(2); plot(h.vchan,rad2bt(h.vchan,p0ALL.rad_allsky)); title('PCRTM allsky');
  figure(3); plot(h.vchan,rad2bt(h.vchan,p0ALL.rad_clrsky)); title('PCRTM clrsky');

  figure(4); plot(h.vchan,rad2bt(h.vchan,p0ALL.rad_clrsky) - rad2bt(h.vchan,p0ALL.rad_allsky)); 
    title('PCRTM : clr - cld')
  figure(5); plot(h.vchan,rad2bt(h.vchan,p0ALL.sarta_clear) - rad2bt(h.vchan,p0ALL.rcalc)); 
    title('SARTA : clr - cld')

  figure(6); plot(h.vchan,rad2bt(h.vchan,p0ALL.rad_clrsky) - rad2bt(h.vchan,p0ALL.sarta_clear)); 
    title('clr : PCRTM - SARTA')
  figure(7); plot(h.vchan,rad2bt(h.vchan,p0ALL.rad_allsky) - rad2bt(h.vchan,p0ALL.rcalc)); 
    title('cld : PCRTM  - SARTA')

  boo = find(h.ichan == 1291);
  figure(8); plot(1:length(p0ALL.stemp),p0ALL.robs1(boo,:),'kx-',...
                  1:length(p0ALL.stemp),p0ALL.sarta_clear(boo,:),'bo-',...
                  1:length(p0ALL.stemp),p0ALL.rcalc(boo,:),'ro-',...
                  1:length(p0ALL.stemp),p0ALL.rad_clrsky(boo,:),'cs-',...
                  1:length(p0ALL.stemp),p0ALL.rad_allsky(boo,:),'ms-')
  hl = legend('obs','Sarta clr','sarta cld','pcrtm clr','pcrtm cld','location','southwest');
  set(hl,'fontsize',10)

end
