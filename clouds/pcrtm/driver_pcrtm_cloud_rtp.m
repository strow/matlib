function [h1ALL,h1aALL,p1ALL,p1aALL] = driver_pcrtm_cloud_rtp(h_inputLVLS,ha,p0ALL_inputLVLS,pa,run_sarta)

% function [h1ALL,h1aALL,p1ALL,p1aALL] = driver_pcrtm_cloud_rtp(h_inputLVLS,ha,p0ALL_inputLVLS,pa,run_sarta)
%     [this is based on xdriver_PCRTM_compute_for_AIRS_spectra_ERAorECMWF.m]
%     [takes about 800 secs for 41 chans/2700 profiles]
%
% takes an input [h_inputLVLS ha p0ALL_inputLVLS pa] which includes cloud structure from (ERA/ECMWF) and
% then runs the PCRTM code; always give estimate of PCRTM clear and PCRTM cloud
%
% p0ALL_inputLVLS is typically 37x12150, or 60x12150, or 90x12150 ... which usually spans surface to 1 mb
%   however, since we will be adding in Xiuhong CO2 profile, and adding in from 1 mb to 0.005 mb, so we
%   will be adding in a few levels!!!
% this means we "pad" p0ALL_inputLVLS to account for this
% 
% can optionally also run SARTA clear and/or SARTA cloud
%
% run_sarta = optional structure argument that says
%   >>> options for PCRTM runs
%     run_sarta.overlap = 1,2,3 for max overlap, random pverlap, max random overlap (suggest 3)
%     run_sarta.ncol0 >= 1 for number of random subcolumns (recommend 50)
%       if run_sarta.ncol0 == -1 then we do ONE column, cloud fraction = 1 == TEST CASE
%       if run_sarta.waterORice = +/-1 we are going to set run_sarta.ncol0 == -1, p.cc = 1 and turn off water or ice clouds
%          so when SARTA is called it turns off appropriate ice or water slab
%          this is test of ONE SLAB CLOUD vs ONE COLUMN CLOUD
%   >>> options for SARTA runs
%     run_sarta.clear           = +/-1 for yes/no, results into prof.sarta_clear
%     run_sarta.cloud           = +/-1 for yes/no, results into prof.sarta_cloudy
%     run_sarta.klayers_code    = string to klayers executable
%     run_sarta.sartaclear_code = string to sarta clear executable
%     run_sarta.sartacloud_code = string to sarta cloud executable
%     run_sarta.randomCpsize        = +1 to randomize BOTH ice (based on Tcld) and water deff
%                                   =  20    (DEFAULT) then water is ALWAYS 20 um (as in PCRTM wrapper), random ice
%                                            (based on Tcld)
%                                   = -1,    then water is MODIS dme, random ice (based on Tcld)
%                                   = +9999, then water is MODIS dme, random ice (based on KN Liou Tcld)
%  run_sarta.iNewVSOrig    = +1;  %% when calling PCRTM_compute_for_AIRS_spectra.m, use new compiled PCRTM_V2.1.exe
%  run_sarta.iUMichCO2     = 0;   %% when calling PCRTM_compute_for_AIRS_spectra.m, use CO2 profile (InputDir/par_constant.dat) and
%                                 %% molindx = 2 and keep scale factor unchanged from what Xiuhong originally gave
%			          %% (1.0135 ==> 385.848*1.0135 = 390.1975)
%     run_sarta.co2ppm          = -1 to use 370 + (yy-2002)*2.2) in pcrtm/sarta
%                               = +x to use user input value     in pcrtm/sarta 
%                               =  0 to use DEFAULT klayers = 385.848 (set in executable by Scott; equivalent to run_sarta.co2ppm = +385.848)
%                   We need to do this as we do not run klayers here, but an internal converter
%                   Note : this value also percolates to klayers/sarta-clear/sarta-cloud, if they are called
% >>> test ONE cloud
%     run_sarta.waterORice = 0 (default, use both clouds)
%                          = +1 turn off ice   clouds, keep water only, set ncol0 = -1, set cc = 1
%                          = -1 turn off water clouds, keep ice   only, set ncol0 = -1, set cc = 1
%       if run_sarta.waterORice = +/-1 we are going to set run_sarta.ncol0 == -1, p.cc = 1 and turn off water or ice clouds
%          so when SARTA is called it turns off appropriate ice or water slab
%          this is test of ONE SLAB CLOUD vs ONE COLUMN CLOUD
%
% Requirements : 
%   p0ALL_inputLVLS must contain ciwc clwc cc from ERA/ECMWF (ie 91xN or 37xN) as well as gas_1 gas_3 ptemp etc
%   the PCRTM code puts in its own particles sizes, cloud fracs and cloud tops based on ciwc,clwc,cc
%   h.ptype = 0 (ie must be levels profile)
%
% Can do ECM (91 levels) or ERA (37 levels)
%
% see below for fixed definitions related to ncol0, overlap
%
%{
quick_test_driver_pcrtm_cloud_rtp
%}
%
%
% testing ONE cloud (driver_pcrtm_cloud_rtp_onecldtest.m)
%   test_onecld_pcrtm
% though remember,
%    simplest way of turning off ice   is set p.ciwc = 0
%    simplest way of turning off water is set p.clwc = 0,
% and then set p.cc = 1

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

if h_inputLVLS.ngas ~= 2
  disp('warning expecting TWO input gases (gas_1 and gas_3) ... if you added in gases 2,4,5 or 6, will replace them with PCRTM profiles ....');
end

if h_inputLVLS.ptype ~= 0
  error('expecting levels profile (h.ptype == 0) ');
end

[ijunk,iajunk,ibjunk] = intersect(h_inputLVLS.glist,1);
if length(ijunk) ~= 1
  error('expecting one of the iunput gases to be gas_1');
end
if h_inputLVLS.glist(iajunk) ~= 1 | h_inputLVLS.gunit(iajunk) ~= 21
  error('expecting gas_1 in g/g');
end

[ijunk,iajunk,ibjunk] = intersect(h_inputLVLS.glist,3);
if length(ijunk) ~= 1
  error('expecting one of the iunput gases to be gas_3');
end
if h_inputLVLS.glist(iajunk) ~= 3 | h_inputLVLS.gunit(iajunk) ~= 21
  error('expecting gas_3 in g/g');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% recommended defaults
ncol0   = 50;  %% number of random overlap clouds

overlap = 1;   %% switch for maximum overlap
overlap = 2;   %% switch for random overlap
overlap = 3;   %% switch for maximum random overlap

if nargin == 4
  run_sarta.waterORice    = 0;   %% keep both clouds
  run_sarta.clear         = -1;
  run_sarta.cloud         = -1;
  %% these are for PCRTM
  run_sarta.ncol0         = 50;
  run_sarta.overlap       = 3;
  run_sarta.randomCpsize = +20;  %% keep Xiangle's ice dme parameterization (based on KN Liou) and 20 um water dme
  run_sarta.co2ppm        = 0;   %% sets default of 385.848 ppm
  run_sarta.iNewVSOrig    = +1;  %% when calling PCRTM_compute_for_AIRS_spectra.m, use new compiled PCRTM_V2.1.exe
  run_sarta.iUMichCO2     = 0;   %% when calling PCRTM_compute_for_AIRS_spectra.m, use CO2 profile (InputDir/par_constant.dat) and
                                 %% molindx = 2 and keep scale factor unchanged from what Xiuhong originally gave
				 %% (1.0135 ==> 385.848*1.0135 = 390.1975)
  addpath ../
  choose_klayers_sarta

elseif nargin == 5
  if ~isfield(run_sarta,'waterORice')
    run_sarta.waterORice    = 0;   %% keep both clouds
  end
  if ~isfield(run_sarta,'co2ppm')
    run_sarta.co2ppm = 0;  %% sets default of 385.848 ppm
  end
  if ~isfield(run_sarta,'clear')
    run_sarta.clear = -1;
  end
  if ~isfield(run_sarta,'cloud')
    run_sarta.cloud = -1;
  end
  %% these are for PCRTM
  if ~isfield(run_sarta,'randomCpsize')
    run_sarta.randomCpsize = +20;
  end
  if ~isfield(run_sarta,'ncol0')
    run_sarta.ncol0 = 50;
  end
  if ~isfield(run_sarta,'overlap')
    run_sarta.overlap = 3;
  end
  if ~isfield(run_sarta,'iNewVSOrig')
    run_sarta.iNewVSOrig = 1;
  end
  if ~isfield(run_sarta,'iUMichCO2')
    run_sarta.iUMichCO2 = 0;
  end

  addpath ../
  choose_klayers_sarta
end

ncol0   = run_sarta.ncol0;
overlap = run_sarta.overlap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mmjunk,nnjunk] = size(p0ALL_inputLVLS.plevs);
fprintf(1,'    >> size of plevs before padding INPUT prof struct = %5i x %5i \n',mmjunk,nnjunk)

[h,p0ALL,p_co2_n2o_co_ch4_pcrtm] = pad_upper_atm(h_inputLVLS,p0ALL_inputLVLS);

[mmjunk,nnjunk] = size(p0ALL.plevs);
fprintf(1,'    >> size of plevs after  padding INPUT prof struct = %5i x %5i \n',mmjunk,nnjunk)

% disp(' >>>>>>>>>>>>>>>>>>> NOT adding in toa stuff in driver_pcrtm_cloud_rtp.m <<<<<<<<<<<<<<<<<<<<<<<<')
% disp(' >>>>>>>>>>>>>>>>>>> NOT adding in toa stuff in driver_pcrtm_cloud_rtp.m <<<<<<<<<<<<<<<<<<<<<<<<')
% h     = h_inputLVLS;
% p0ALL = p0ALL_inputLVLS;

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

if run_sarta.waterORice == -1
  %% removes water, keeps ice clouds
  disp(' >>> WARNING remove water, keep ice clouds in PCRTM <<<<<<<<');
  disp(' >>> WARNING remove water, keep ice clouds in PCRTM <<<<<<<<');  
  p0ALL.clwc = 0 * p0ALL.clwc;
  run_sarta.ncol0 = -1;
  ncol0 = -1;
  p0ALL.cc = ones(size(p0ALL.cc));  
elseif run_sarta.waterORice == +1
  %% keeps water, removes ice clouds
  disp(' >>> WARNING remove ice, keep water clouds in PCRTM <<<<<<<<');
  disp(' >>> WARNING remove ice, keep water clouds in PCRTM <<<<<<<<');    
  p0ALL.ciwc = 0 * p0ALL.ciwc;
  run_sarta.ncol0 = -1;
  ncol0 = -1;
  p0ALL.cc = ones(size(p0ALL.cc));  
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

iChunk = 100;  %% speed up the code  by breaking input profiles into chunks, don't change this (code really slows down!)

clear p

p1ALL  = p0ALL;
p0ALLX = p0ALL;

% 12784 * 86400 + 27 = 1.1045e+09;
if nanmean(p0ALL.rtime) > 1e9
  %% /asl/matlab2012/airs/readers/xreadl1b_all.m
  [yy,mm,dd,hh] = tai2utc(p0ALL.rtime - (12784 * 86400 + 27));
else
  [yy,mm,dd,hh] = tai2utc(p0ALL.rtime);
end

if run_sarta.randomCpsize == +9999 |  run_sarta.randomCpsize == -1
  modis_waterDME = modisL3_map_rtp_cloudprofile(p0ALLX);
else
  modis_waterDME = [];
  modis_waterDME = -9999 * ones(size(p0ALL.stemp));
end

iIndMax = ceil(length(p0ALL.xtrack)/iChunk);

fprintf(1,'num of input profiles = %4i will be processed in %4i chunks \n',length(p0ALL.xtrack),iIndMax);

for iInd = 1 : iIndMax

  inds = (1:iChunk) + (iInd-1)*iChunk;
  inds = intersect(1:length(p0ALL.xtrack),inds);

  %p0 = index_subset(inds,p0ALLX); 
  %[h,ha,p,pa] = rtpgrow(h,ha,p0,pa);
  %[h,p] = subset_rtp_clouds(h,p0ALLX,[],[],inds);  
  [h,p] = subset_rtp_allcloudfields(h,p0ALLX,[],[],inds);

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

  %%% zen_ang = double(p.scanang);      %% orig, gives very large average clr sky biases between SARTA and PCRTM
  zen_ang = double(abs(p.satzen));  %% new and agrees much better with SARTA clear sky, Sergio 08/19/2015

  if run_sarta.co2ppm == -1
    %co2     = ones(size(p.stemp)) .* (370 + (yy(inds')-2002)*2.2);  
    deltaT = (yy(inds')-2002) + (mm(inds')-1)/12 + dd(inds')/30/12;
    co2    = ones(size(p.stemp)) .* (370 + deltaT*2.2);
  elseif run_sarta.co2ppm == 0
    co2    = ones(size(p.stemp)) * 385.848;   %% and pcrtm essentially does not touch anything in PCRTM_compute_for_AIRS_spectra
  elseif run_sarta.co2ppm > 0
    co2    = ones(size(p.stemp)) * run_sarta.co2ppm;
  end
  co2_all(inds) = co2;
  
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

  iDoCalcPCRTM = -1;  %% for fast SARTA debugging
  iDoCalcPCRTM = +1;  %% what we want : PCRTM calcs!!!
  
  if iDoCalcPCRTM > 0
    %% note that internally this soubroutine uses abs(ncol) so if we use ncol0 = -1, we have ONE column  
    [rad_allsky rad_clrsky tmpjunk rad_allsky_std sarta_gas_2_6] = ...
                                    PCRTM_compute_for_AIRS_spectra(nboxes,nlev, ncol, overlap, ...
                                                           P, WCT, ICT, cc, TT, q, o3, Ps, Ts, ...
                                                           sfctype, efreq, emis, ...
                                                           zen_ang, co2, parname, ppath, ...
                                                           run_sarta.randomCpsize, modis_waterDME(inds),run_sarta);
%    [rad_allsky2 rad_clrsky tmpjunk2] = PCRTM_compute_for_AIRS_spectra_V2(nboxes,nlev, ncol, overlap, ...
%                                                           P, WCT, ICT, cc, TT, q, o3, Ps, Ts, ...
%                                                           sfctype, efreq,emis, ...
%                                                           zen_ang, co2, parname,ppath,run_sarta);
  else
    disp('skipped PCRTM/MRO ... going straight on to SARTA 2S')
    sarta_gas_2_6.co2 = 385.848;
    sarta_gas_2_6.ch4 = 1.843;
    
    rad_allsky     = zeros(2378,length(p.stemp));
    rad_clrsky     = zeros(2378,length(p.stemp));
    rad_allsky_std = zeros(2378,length(p.stemp));
    
    tmpjunk.totalODice  = zeros(1,length(p.stemp));
    tmpjunk.meanDMEice  = zeros(1,length(p.stemp));    
    tmpjunk.maxCTOPice  = zeros(1,length(p.stemp));
    tmpjunk.totalODiceX = zeros(1,length(p.stemp));
    tmpjunk.lvlODice    = zeros(max(p.nlevs),length(p.stemp));
    
    tmpjunk.totalODwater  = zeros(1,length(p.stemp));
    tmpjunk.meanDMEwater  = zeros(1,length(p.stemp));    
    tmpjunk.maxCTOPwater  = zeros(1,length(p.stemp));
    tmpjunk.totalODwaterX = zeros(1,length(p.stemp));
    tmpjunk.lvlODwater    = zeros(max(p.nlevs),length(p.stemp));    
  end

  %% we know we are adding on a few level from the 101 level CO2 profile, into the input ECM/ERA/MERRA profile
  %% be be a little careful
  p_orig_levels = p;

  if run_sarta.clear > 0
    get_sarta_clear2
  end

  if run_sarta.cloud > 0
    %% probably will re-do sarta clear calcs but checks there is zero difference with above clear calcs
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/rtptools
h1ALL       = h;
h1ALL.ngas  = 6;
h1ALL.glist = [ 1  2  3  4  5  6]';
h1ALL.gunit = [21 10 21 10 10 10]';
h1ALL.pmin  = min(min(p1ALL.plevs));
h1aALL      = ha;

% now overwrite p.rcalc and replace with pcrtm calcs
p1ALL.rcalc       = p1ALL.rad_allsky;
p1ALL.rcalc_std   = p1ALL.rad_allsky_std;

%% in ppmv, GUNITS = 10
for ij = 1 : length(p1ALL.stemp)
  p1ALL.gas_2(:,ij) = p_co2_n2o_co_ch4_pcrtm.gas_2(:,ij) * co2_all(ij)/385.848;
  p1ALL.gas_4(:,ij) = p_co2_n2o_co_ch4_pcrtm.gas_4(:,ij)
  p1ALL.gas_5(:,ij) = p_co2_n2o_co_ch4_pcrtm.gas_5(:,ij)    
  p1ALL.gas_6(:,ij) = p_co2_n2o_co_ch4_pcrtm.gas_6(:,ij) * sarta_gas_2_6.ch4/1.843;
end

p1aALL = pa;
p1aALL = set_attr(p1aALL,'sarta_clear',  run_sarta.sartaclear_code);
p1aALL = set_attr(p1aALL,'sarta_cloud',  run_sarta.sartacloud_code);
p1aALL = set_attr(p1aALL,'iUMichCO2',    num2str(run_sarta.iUMichCO2));
p1aALL = set_attr(p1aALL,'cumsum',       num2str(run_sarta.cumsum));
p1aALL = set_attr(p1aALL,'ncol0',        num2str(run_sarta.ncol0));
p1aALL = set_attr(p1aALL,'co2ppm',       num2str(run_sarta.co2ppm));
p1aALL = set_attr(p1aALL,'ForceNewSlabs',num2str(run_sarta.ForceNewSlabs));
