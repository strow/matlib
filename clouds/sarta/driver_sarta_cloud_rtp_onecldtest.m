function [prof,index] = driver_sarta_cloud_rtp_onecldtest(h,ha,p,pa,run_sarta,waterORice)

%% modelled on MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE_CFRAC/readecmwf91_nearest_gasNcloud_slabprof.m
%% also see /asl/rtp_prod/airs/rtp/create_rcalc_ecm_cld_allfov.m
%
% function prof = driver_sarta_cloud_rtp_onecldtest(h,ha,p,pa,run_sarta,waterORice)
% takes an input [h ha p pa] which incudes cloud structure from (ERA/ECMWF) and
% then runs the SARTA code. BUT it also dumps out selected profiles, saving only those which are WATER or ICE CLD only
% depending on setting of waterORice
%
% run_sarta = optional structure argument that says
%   >>> options for SARTA runs
%     run_sarta.clear = +/-1 for yes/no, results into prof.clearcalc
%     run_sarta.cloud = +/-1 for yes/no, results into prof.rcalc
%     run_sarta.cumsum = 0 -- 1 to set cloud pressure based on cumulative sum, -ve for just go with "ecmwf2sarta" results
%     run_sarta.klayers_code = string to klayers
%     run_sarta.sartaclear_code = string to sarta clear executable
%     run_sarta.sartacloud_code = string to sarta cloud executable
%     run_sarta.ice_water_separator = set all ciwc/clwc to ice above this, water below this 
%        (default = -1, use ciwc/clwc structures as is)
%     run_sarta.randomCpsize        = +1 or 0 to turn on/off randomizing of ice and water deff
%                                      if 0, then water is ALWAYS 20 um (as in PCRTM wrapper)
%
% Requirements : 
%   p must contain ciwc clwc cc from ERA/ECMWF (ie 91xN or 37xN) as well as gas_1 gas_3 ptemp etc
%   this code puts in its own particles sizes, cloud fracs and cloud tops based on ciwc,clwc,cc
%   h.ptype = 0 (ie must be levels profile)
%
% Can do ECM (91 levels) or ERA (37 levels)
%
% testing : 
%   test_onecld_sarta
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
%% defaults

run_sarta.cfrac = 1;

if nargin == 4
  run_sarta.clear = -1;
  run_sarta.cloud = +1;
  run_sarta.cumsum = -1;
  %run_sarta.klayers_code = '/asl/packages/klayers/Bin/klayers_airs';
  %run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  run_sarta.ice_water_separator = -1;
  run_sarta.randomCpsize        = +1;  %% keep randomizing dme for ice and water

  waterORice = +1; % keep only water clds
elseif nargin == 5 | nargin == 6
  waterORice = +1; % keep only water clds
  if ~isfield(run_sarta,'randomCpsize')
    run_sarta.randomCpsize = +1;
  end
  if ~isfield(run_sarta,'ice_water_separator')
    run_sarta.ice_water_separator = -1;
  end
  if ~isfield(run_sarta,'clear')
    run_sarta.clear = -1;
  end
  if ~isfield(run_sarta,'cloud')
    run_sarta.cloud = -1;
   end
  if ~isfield(run_sarta,'cumsum')
    run_sarta.cumsum = -1;
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

% Min allowed cloud fraction
cmin = 0.0001;

% Max allowed cngwat[1,2]
cngwat_max = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(p,'ciwc')
  error('driver_pcrtm_cloud_rtp.m requires ciwc');
elseif ~isfield(p,'clwc')
  error('driver_pcrtm_cloud_rtp.m requires clwc');
elseif ~isfield(p,'cc')
  error('driver_pcrtm_cloud_rtp.m requires cc');
elseif h.ptype ~= 0
  error('driver_pcrtm_cloud_rtp.m requires LEVELS profiles (h.ptype = 0)');
end

if ~isfield(p,'cfrac')
  %% need random cfracs
  disp('>>>>>>>> warning : need random cfracs .... initializing')
  %% want to make sure there are NO zeros cfrac
  p.cfrac = 0.50*(rand(size(p.stemp)) + rand(size(p.stemp))) ;
end

if run_sarta.ice_water_separator > 0
  disp('>>>>>>>> warning : setting SEPRATOR for ice and water .... initializing')
  p = convert_ice_water_separator(p,run_sarta.ice_water_separator);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
head  = h;

nlev     = ceil(mean(p.nlevs));
nlev_std = (std(double(p.nlevs)));

if nlev_std > 1e-3
  error('oops : code assumes ERA (37 levs) or ECMWF (91 levs) or other constant numlevs model')
end

if h.ptype ~= 0
  error('need levels input!')
end

load airsheights.dat
load airslevels.dat

tic
[prof,profX] = ecmwfcld2sartacld(p,nlev,run_sarta.cumsum,airslevels,airsheights);   
                                            %% figure the two slab cloud profile info here, using profX
                                            %% this then puts the info into "prof" by calling put_into_prof w/in routine

prof = put_into_V201cld_fields(prof);    %% puts cloud info from above into rtpv201 fields 
  prof.ctype  = double(prof.ctype);
  prof.ctype2 = double(prof.ctype2);

%% sets fracs and particle effective sizes eg cfrac2
prof = set_fracs_deffs(head,prof,profX,cmin,cngwat_max,run_sarta.cfrac,run_sarta.randomCpsize);

if run_sarta.cumsum > 0
  prof = reset_cprtop(prof);
end

prof = check_for_errors(prof,run_sarta.cfrac);           %% see if there are possible pitfalls

clear profX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now turn off ice, change water!!!!!!!!
[prof,index] = only_waterORice_cloud(h,prof,waterORice);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  tic
  get_sarta_clear;
  toc
  prof.rclearcalc = profRX2.rcalc;
end
if run_sarta.cloud > 0 
  disp('running SARTA cloud')
  tic
  get_sarta_cloud;
  toc
  prof.rcalc = profRX2.rcalc;
else
  disp('you did not ask for SARTA cloudy to be run; not changing p.rcalc')
end

tnow = toc;
fprintf(1,'TOTAL : %8.6f minutes to process \n',tnow/60);

 
