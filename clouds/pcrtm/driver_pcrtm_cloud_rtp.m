function [p1ALL] = driver_pcrtm_cloud_rtp(h,ha,p0ALL,pa,run_sarta)

%     [this is based on xdriver_PCRTM_compute_for_AIRS_spectra_ERAorECMWF.m]
%     [takes about 800 secs for 41 chans/2700 profiles]
%
% takes an input [h ha p0ALL pa] which incudes cloud structure from (ERA/ECMWF) and
% then runs the PCRTM code; always give estimate of PCRTM clear and PCRTM cloud
%
% can optionally also run SARTA clear and/or SARTA cloud
%
% run_sarta = optional structure argument that says
%   >>> options for PCRTM runs
%     run_sarta.overlap = 1,2,3 for max overlap, random pverlap, max random overlap (suggest 3)
%     run_sarta.ncol0 >= 1 for number of random subcolumns (recommend 50)
%       if run_sarta.ncol0 == -1 then we do ONE column, cloud fraction = 1 == TEST CASE
%   >>> options for SARTA runs
%     run_sarta.clear = +/-1 for yes/no, results into prof.sarta_clear
%     run_sarta.cloud = +/-1 for yes/no, results into prof.sarta_cloudy
%     run_sarta.klayers_code = string to klayers executable
%     run_sarta.sartaclear_code = string to sarta clear executable
%     run_sarta.sartacloud_code = string to sarta cloud executable
%
% Requirements : 
%   p0ALL must contain ciwc clwc cc from ERA/ECMWF (ie 91xN or 37xN) as well as gas_1 gas_3 ptemp etc
%   the PCRTM code puts in its own particles sizes, cloud fracs and cloud tops based on ciwc,clwc,cc
%   h.ptype = 0 (ie must be levels profile)
%
% Can do ECM (91 levels) or ERA (37 levels)
%
% see below for fixed definitions related to ncol0, overlap
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fileparts(mfilename('fullpath')); % current directory
base_dir1 = fileparts(base_dir);  % dir:  ../
base_dir2 = fileparts(base_dir1); % dir:  ../../

% Airs Matlab utility box
addpath([base_dir2 '/gribtools'])
addpath([base_dir2 '/rtptools'])
addpath([base_dir2 '/aslutil'])
addpath([base_dir2 '/science'])
addpath([base_dir2 '/h4tools'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
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
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% recommended defaults
ncol0   = 50;  %% number of random overlap clouds

overlap = 1;   %% switch for maximum overlap
overlap = 2;   %% switch for random overlap
overlap = 3;   %% switch for maximum random overlap

if nargin == 4
  run_sarta.clear = -1;
  run_sarta.cloud = -1;
  run_sarta.ncol0 = 50;
  run_sarta.overlap = 3;
elseif nargin == 5
  if ~isfield(run_sarta,'clear')
    run_sarta.clear = -1;
  end
  if ~isfield(run_sarta,'cloud')
    run_sarta.cloud = -1;
  end
  if ~isfield(run_sarta,'overlap')
    run_sarta.overlap = 3;
  end
  if ~isfield(run_sarta,'ncol0')
    run_sarta.ncol0 = 50;
  end
  if ~isfield(run_sarta,'klayers_code')
    %run_sarta.klayers_code = '/asl/packages/klayers/Bin/klayers_airs'; 
    run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  end   
  if ~isfield(run_sarta,'sartaclear_code')
    run_sarta.sartaclear_code = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';
  end
  if ~isfield(run_sarta,'sartacloud_code')
    %run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
    run_sarta.sartacloud_code = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  end
end

ncol0 = run_sarta.ncol0;
overlap = run_sarta.overlap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(p0ALL,'ciwc')
  error('driver_pcrtm_cloud_rtp.m requires ciwc');
elseif ~isfield(p0ALL,'clwc')
  error('driver_pcrtm_cloud_rtp.m requires clwc');
elseif ~isfield(p0ALL,'cc')
  error('driver_pcrtm_cloud_rtp.m requires cc');
elseif h.ptype ~= 0
  error('driver_pcrtm_cloud_rtp.m requires LEVELS profiles (h.ptype = 0)');
end

iChunk = 100;  %% speed up the code  by breaking input profiles into chunks, don't change this (code really slows down!)

clear p

p1ALL = p0ALL;
p0ALLX = p0ALL;

[yy,mm,dd,hh] = tai2utc(p0ALL.rtime);

iIndMax = ceil(length(p0ALL.xtrack)/iChunk);

fprintf(1,'num of input profiles = %4i will be processed in %4i chunks \n',length(p0ALL.xtrack),iIndMax)

for iInd = 1 : iIndMax

  inds = (1:iChunk) + (iInd-1)*iChunk;
  inds = intersect(1:length(p0ALL.xtrack),inds);

  %p0 = index_subset(inds,p0ALLX); 
  %[h,ha,p,pa] = rtpgrow(h,ha,p0,pa);
  [h,p] = subset_rtp_clouds(h,p0ALLX,[],[],inds);

  nboxes = length(p.stemp);  
  [nlev,nprof] = size(p.clwc);
  ncol = ncol0;

  P   = double(p.plevs);  %% 1 mb = 1 hPa
  
  WCT = double(p.clwc);   %% cloud liquid water content in kg/kg
  ICT = double(p.ciwc);   %% cloud ice    water content in kg/kg
  cc  = double(p.cc);     %% cloud fraction
  if ncol0 == -1
    disp('FORCE CFRAC = 1 at all levels : TEST CASE');
    yes_cld = find(cc > eps);
    cc(yes_cld) = 1;
  end

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
  co2     = ones(size(p.stemp)) .* (370 + (yy(inds')-2002)*2.2);

  % use_Xiuhong  %% for debug default 2012/05/01  00:00-01:00 UTC

  ix = inds(1);
  yymmddhhstr = [num2str(yy(ix)) '.' num2str(mm(ix),'%02d') '.' num2str(dd(ix),'%02d') '.' num2str(hh(ix))];

  %parname = ['/strowdata1/s1/sergio/PCRTM_XIANGLEI/CLEANCOPY/JUNK/pcrtm' yymmddhhstr '.tmp'];
  parname = mktemp(['pcrtm' yymmddhhstr]);

  ppath = '/strowdata1/s1/sergio/PCRTM_XIANGLEI/NEWVERS/PCRTM_V2.1_for_AIRS/code_changed/Run/';
  ppath = '/asl/s1/sergio/PCRTM_XIANGLEI/NEWVERS/PCRTM_V2.1_for_AIRS/code_changed/Run/';
  parnameout = [parname '.out'];

  %whos P WCT ICE cc TT q o3 Ps Ts sfctype efreq emis zen_ang co2 
  fprintf(1,'making PCRTM input file %s for iChunk %3i of %3i \n',parname,iInd,iIndMax)

  [rad_allsky rad_clrsky tmpjunk] = PCRTM_compute_for_AIRS_spectra(nboxes,nlev, ncol, overlap, ...
                                                           P, WCT, ICT, cc, TT, q, o3, Ps, Ts, ...
                                                           sfctype,efreq,emis, ...
                                                           zen_ang,co2, parname,ppath);
%    [rad_allsky2 rad_clrsky tmpjunk2] = PCRTM_compute_for_AIRS_spectra_V2(nboxes,nlev, ncol, overlap, ...
%                                                           P, WCT, ICT, cc, TT, q, o3, Ps, Ts, ...
%                                                           sfctype,efreq,emis, ...
%                                                           zen_ang,co2, parname,ppath);

  if run_sarta.clear > 0
    get_sarta_clear2
  end

  if run_sarta.cloud > 0
    get_sarta_cloud2
  end

  get_extra_p1_params_inds

  rmer = ['!/bin/rm ' parname ' ' parnameout];
  eval(rmer);

end

if run_sarta.cloud > 0
  disp('added on sarta cloud calcs');
end
if run_sarta.clear > 0
  disp('added on sarta clear calcs');
end

% now overwrite p.rcalc and replace with pcrtm calcs
p1ALL.rcalc = p1ALL.rad_allsky;

