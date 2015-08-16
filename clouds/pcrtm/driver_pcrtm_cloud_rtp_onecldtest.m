function [p1ALL] = driver_pcrtm_cloud_rtp_onecldtest(h,ha,p0ALL,pa,run_sarta,waterORice)

%     [this is based on xdriver_PCRTM_compute_for_AIRS_spectra_ERAorECMWF.m]
%     [takes about 800 secs for 41 chans/2700 profiles]
%
% function prof = driver_pcrtm_cloud_rtp_onecldtest(h,ha,p,pa,run_sarta,waterORice)
% takes an input [h ha p0ALL pa] which incudes cloud structure from (ERA/ECMWF) and
% then runs the PCRTM code; always give estimate of PCRTM clear and PCRTM cloud.
% BUT it also dumps out selected profiles, saving only those which are WATER or ICE CLD only
% depending on setting of waterORice
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
% testing
%   test_onecld_pcrtm
% though remember,
%    simplest way of turing off ice   is set p.ciwc = 0
%    simplest way of turing off water is set p.clwc = 0,
% and then set p.cc = 1
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
%% recommended defaults
ncol0   = 50;  %% number of random overlap clouds

overlap = 1;   %% switch for maximum overlap
overlap = 2;   %% switch for random overlap
overlap = 3;   %% switch for maximum random overlap

if nargin == 4
  run_sarta.clear = -1;
  run_sarta.cloud = -1;
  run_sarta.ncol0         = 50;
  run_sarta.overlap       = 3;
  run_sarta.randomCpsize = +20;  %% keep Xiangle's ice dme parmerization (based on KN Liou) and 20 um water dme
  run_sarta.co2_ppm       = 0;   %% sets default of 385.848 ppm  
  addpath ../
  choose_klayers_sarta

elseif nargin == 5
  if ~isfield(run_sarta,'co2_ppm')
    run_sarta.co2_ppm = 0;  %% sets default of 385.848 ppm  
  end
  if ~isfield(run_sarta,'clear')
    run_sarta.clear = -1;
  end
  if ~isfield(run_sarta,'cloud')
    run_sarta.cloud = -1;
  end
  if ~isfield(run_sarta,'overlap')
    run_sarta.overlap = 3;
  end
  if ~isfield(run_sarta,'randomCpsize')
    run_sarta.randomCpsize = +20;
  end  
  if ~isfield(run_sarta,'ncol0')
    run_sarta.ncol0 = 50;
  end

  addpath ../
  choose_klayers_sarta
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

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if waterORice == +1
  p0ALL.clwc = 0 * p0ALL.clwc;
elseif waterORice == -1
  p0ALL.ciwc = 0 * p0ALL.ciwc;
end
run_sarta.ncol0 = -1;
ncol0 = -1;
p0ALL.cc = ones(size(p0ALL.cc));
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

iChunk = 100;  %% speed up the code  by breaking input profiles into chunks, don't change this (code really slows down!)

clear p

p1ALL = p0ALL;
p0ALLX = p0ALL;

% 12784 * 86400 + 27 = 1.1045e+09;
if nanmean(p0ALL.rtime) > 1e9
  %% /asl/matlab2012/airs/readers/xreadl1b_all.m
  [yy,mm,dd,hh] = tai2utc(p0ALL.rtime - (12784 * 86400 + 27));
else
  [yy,mm,dd,hh] = tai2utc(p0ALL.rtime);
end

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

  if run_sarta.co2_ppm == -1
    %co2     = ones(size(p.stemp)) .* (370 + (yy(inds')-2002)*2.2);
    deltaT = (yy(inds')-2002) + (mm(inds')-1)/12 + dd(inds')/30/12;
    co2    = ones(size(p.stemp)) .* (370 + deltaT*2.2);    
  elseif run_sarta.co2_ppm == 0
    co2     = ones(size(p.stemp)) * 385.848;
  elseif run_sarta.co2_ppm > 0
    co2     = ones(size(p.stemp)) * run_sarta.co2_ppm;
  end
  co2_all(inds) = co2;

  % use_Xiuhong  %% for debug default 2012/05/01  00:00-01:00 UTC

  ix = inds(1);
  yymmddhhstr = [num2str(yy(ix)) '.' num2str(mm(ix),'%02d') '.' num2str(dd(ix),'%02d') '.' num2str(hh(ix))];

  %parname = ['/strowdata1/s1/sergio/PCRTM_XIANGLEI/CLEANCOPY/JUNK/pcrtm' yymmddhhstr '.tmp'];
  parname = mktemp(['pcrtm' yymmddhhstr]);

  ppath = '/strowdata1/s1/sergio/PCRTM_XIANGLEI/NEWVERS/PCRTM_V2.1_for_AIRS/code_changed/Run/';
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

