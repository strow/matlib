function [prof,index] = driver_sarta_cloud_rtp_onecldtest(h,ha,p,pa,run_sarta)

%% code can esily be duplicated using  driver_sarta_cloud_rtpm
%%
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
%     run_sarta.co2ppm          = -1 to use 370 + (yy-2002)*2.2) in pcrtm/sarta
%                               = +x to use user input value     in pcrtm/sarta 
%                               =  0 to use DEFAULT klayers = 385 (set in executable by Scott; equivalent to run_sarta.co2ppm = +385)
%                   this is done to keep it consistent with PCRTM
%                   however, also have to make sure this is only enabled if h.glist does NOT include gas_2
%      run_sarta.ForceNewSlabs  = -1 (default) to keep any slab clouds that are input, as they are
%                               = +1           to force new slabs to be derived from clwc,ciwc,cc
%      run_sarta.waterORice     = +1/-1        to do water/ice clouds
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
check_sarta_cloud_rtp_defaults
run_sarta.cfrac = 1;  %% set cfrac == 1

%% turn profiles into slabs
main_code_to_make_slabs

%% add on co2
prof_add_co2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now turn off ice, change water!!!!!!!!
[prof,index] = only_waterORice_cloud(h,prof,run_sarta.waterORice);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  tic
  get_sarta_clear;
  toc
  prof.rclearcalc = profRX2.rcalc;
else
  disp('you did not ask for SARTA clear to be run; not changing p.sarta_rclearcalc')
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

 
