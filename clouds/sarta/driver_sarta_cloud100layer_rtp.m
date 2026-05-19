function prof = driver_sarta_cloud100layer_rtp(h,ha,p,pa,run_sarta)

%% modelled on MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE_CFRAC/readecmwf91_nearest_gasNcloud_slabprof.m
%% also see /asl/rtp_prod/airs/rtp/create_rcalc_ecm_cld_allfov.m
%
% function prof = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta)
% takes an input [h ha p pa] which incudes cloud structure from (ERA/ECMWF) and
% then runs the SARTA code
%
% run_sarta = optional structure argument that says
% >>> options for SARTA runs
%     run_sarta.clear = +/-1 for yes/no, results into prof.clearcalc
%     run_sarta.cloud = +/-1 for yes/no, results into prof.rcalc
%     >>> no need for run_sarta.cumsum <<<<<
%     run_sarta.ncol >= 1              : number of subcolumns (if cfrac == 1)
%     run_sarta.overlap = 1,2,3        : type of cloud overlap if ncol > 0
%     run_sarta.cfrac < 0              : use random
%                     > 0 to < 1       : use fixed amount specified by user
%     run_sarta.klayers_code        = string to klayers
%     run_sarta.sartaclear_code     = string to sarta clear executable
%     run_sarta.sartacloud_code     = string to sarta cloud executable
%     run_sarta.ice_water_separator = set all ciwc/clwc to ice above this, water below this
%                                   = -1;  %% DEFAULT = -1, DO NOT CALL convert_ice_water_separator, used since 2010(???) onwards
%                                                           cloud_combine_main_code uses ISCCP ie ice above 440 mb, water below 440 mb
%                                   = 0;   %% do not separate out ciwc and clwc by pressure; ie believe the NWP are correct
%                                                           do NOT CALL convert_ice_water_separator
%                                                           cloud_combine_main_code uses nothing
%                                   = +1;  %% use quadratic X = [-60 0 +60]; Y = [6 9 6]; P = polyfit(X,Y,2); X1=[p.rlat];Y1=polyval(P,X1); according to IPCC AR5
%                                                           DO NOT CALL convert_ice_water_separator
%                                                           cloud_combine_main_code uses ISCCP ie ice above Y1 mb, water below Y1 mb
%                                   >>>>>>> these first alter the ciwc and clwc profiles <<<<<<<<<
%                                   = +2;  %% use quadratic X = [-60 0 +60]; Y = [6 9 6]; P = polyfit(X,Y,2); X1=[p.rlat];Y1=polyval(P,X1); according to IPCC AR5
%                                                           DO CALL convert_ice_water_separator
%                                                           cloud_combine_main_code uses ISCCP ie ice above Y1 mb, water below Y1 mb
%                                   = +440;%% or similar [100 -- 1000] ... almost same as DEFAULT = -1, but
%                                                           DO CALL convert_ice_water_separator
%                                                           cloud_combine_main_code uses ISCCP ie ice above 440 mb, water below X mb
%     run_sarta.randomCpsize        = +1 (default) to randomize BOTH ice (based on Tcld) and water deff
%                                      20,   then water is ALWAYS 20 um (as in PCRTM wrapper), random ice
%                                      (based on Tcld)
%                                      -1, then water is MODIS dme, random ice (based on Scott Tcld)
%                                      9999, then water is MODIS dme, random ice (based on KNLiou Tcld)
%     run_sarta.co2ppm          = -1 to use 370 + (yy-2002)*2.2) in pcrtm/sarta
%                               = +x to use user input value     in pcrtm/sarta 
%                               =  0 to use DEFAULT klayers = 385 (set in executable by Scott; equivalent to run_sarta.co2ppm = +385)
%                   this is done to keep it consistent with PCRTM
%                   however, also have to make sure this is only enabled if h.glist does NOT include gas_2
%%
% Requirements : 
%   p must contain ciwc clwc cc from ERA/ECMWF (ie 91xN or 37xN) as well as gas_1 gas_3 
%   ptemp etc
%   this code puts in its own particles sizes, cloud fracs and cloud tops based on ciwc,clwc,cc
%   h.ptype = 0 (ie must be levels profile)
%
% Can do arbitrary levels eg ECM (91 levels) or ERA (37 levels)

% can do "driver_sarta_cloud_rtp_onecldtest", but remember
%    simplest way of turing off ice   is set p.ciwc = 0
%    simplest way of turing off water is set p.clwc = 0,
% and then set p.cc = 1

%
% Written by Sergio DeSouza-Machado (with a lot of random cloud frac and dme by Scott Hannon)
%
% BASICALLY SAME AS driver_sarta_cloud_rtp.m EXCEPT it does NOT need run_sarta.cumsum
% and instead needs run_sarta.ncol >= 1               : number of subcolumns (if cfrac == 1)
%
% updates
%  04/27/2014 : introduced run_sarta.ice_water_separator + 1 option, which is
%               quadratic as fcn of latitude separator, acording to IPCC AR5 report (Ch 7, Fig 7.5)
%               X = [-60 0 +60]; Y = [6 9 6]; P = polyfit(X,Y,2); X1=[-90:5:+90];Y1=polyval(P,X1);
%               above certain pressure = ice      below certain pressure = water
%  08/18/2013 : introduced run_sarta.randomCpsize, default = +1 to keep randomizing deff; 
%                                                             20 to keep ice = ice(T,pcrtm), water = 20 um
%                                                             -1 to use MODIS water DME
%  08/17/2013 : if cfrac does not exist, fix random initialization so as to make sure cfrac > 0 always
%  08/17/2013 : making more extensive use of run_sarta.cfrac by introducing it in set_fracs_deffs.m and check_for_errors.m
%  08/15/2013 : introduced run_sarta.ice_water_separator, so that everything above certain pressure == ice; 
%               below certain pressure = water
%  04/22/2013 : introduced run_sarta.cfrac so that we can eg have cfrac = 1
%  04/04/2013 : introduced/tweaked do_the_reset_cprtop_cloudOD.m
%  04/04/2013 : mega updates in check_for_errors.m
%
%{
% testing the code for 100layer clouds, random overlap
  [h,ha,p,pa] = rtpread('/asl/data/rtprod_airs/2012/05/01/cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp');
  [h,ha,p,pa] = rtpgrow(h,ha,p,pa);
  run_sarta.clear = +1;
  run_sarta.cloud = +1;
  run_sarta.ncol  = 50;
  run_sarta.cfrac = 1;
  run_sarta.ice_water_separator = +1;    %% model 1, use X = [-60 0 +60]; Y = [6 9 6];
                                         %% P = polyfit(X,Y,2); X1=[-90:5:+90];Y1=polyval(P,X1); from IPCC AR5 report (Chapter 7, Fig.7.5):
					 %% see https://www.gfdl.noaa.gov/clouds-climate-initiative/
  run_sarta.ice_water_separator = 440;   %% this is from ISCCP (ice clouds above 440 mb)
  tic
  p1 = driver_sarta_cloud100layer_rtp(h,ha,p,pa,run_sarta);
  toc
  rtpwrite('/asl/data/rtprod_airs/2012/05/01/pcrtm_cld_ecm_41ch.airs_ctr.2012.05.01.10_sarta.rtp',h,ha,p1,pa);
%}
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
narginx = nargin;
if nargin == 4
  run_sarta = struct;
end  
[p,run_sarta,otherstuff] = check_sarta_cloud_rtp_defaults(run_sarta,h,p,narginx);

cmin = otherstuff.cmin;             %% min allowed cfrac
cngwat_max = otherstuff.cngwat_max; %% max allowed cngwat
iDebugMain = otherstuff.iDebugMain; %% to debug or not???

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do 2 slabclouds

%% turn profiles into slabs, waste of time
%% prof = main_code_to_make_slabs(h,ha,p,pa,run_sarta,iDebugMain,otherstuff);
prof = p;

%% add on co2
prof = prof_add_co2(h,prof,run_sarta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(prof,'rcalc')
  prof = rmfield(prof,'rcalc');
end
if isfield(prof,'rcld')
  prof = rmfield(prof,'rcld');
end

run_sarta.Slab_or_100layer = -1;     %% now explicitly set this
if run_sarta.Slab_or_100layer == -1  %% which it WILL be, given line above
  if ~isfield(run_sarta,'ncol')
    run_sarta.ncol    =  1;  %% number of columns for 100 layer cloud code
    run_sarta.ncol    = 25;  %% number of columns for 100 layer cloud code
  end
  if ~isfield(run_sarta,'overlap')
    run_sarta.overlap = +3;  %% maximal random overlap for 100 layer cloud code
  end
  %% see /home/sergio/klayersV205/Src_rtpV201_100layercloudamountsize
  run_sarta.klayers_code = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';
  %% see /home/sergio/SARTA_CLOUDY/v108_Src_rtpV201_pclsam_slabcloud_hg3_100layerNEW
  run_sarta.sartacloud_code = ...
      '/home/sergio/SARTA_CLOUDY/BinV201/WORKS_Dec2015/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';
end

%% now essentially IGNORE the above and instead do the cloud PROFILES stuff
[h,prof] = reset_cloud_slab_with_cloud_profile(h,prof,run_sarta.cfrac);

%% and add on aux info, such as OD etc
load airsheights.dat
load airslevels.dat

airsheights = flipud(airsheights);
airslevels  = flipud(airslevels);
playsN = airslevels(1:end-1)-airslevels(2:end);
playsD = log(airslevels(1:end-1)./airslevels(2:end));
airslayers = playsN./playsD;

prof = cloudOD_for100layer(prof,run_sarta.cumsum/100,airslevels,airslayers,airsheights);

%%%%%%%%%%%%%%%%%%%%%%%%%

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  prof = get_sarta_clear(h,ha,prof,pa,run_sarta);
else
  disp('you did not ask for SARTA clear to be run; not changing p.rcalc')
end

if run_sarta.cloud > 0 
  if run_sarta.ncol == 1
    %% all you have to do is run ONCE
    disp('running SARTA cloud ONCE')
    prof = get_sarta_cloud100layer(h,ha,prof,pa,run_sarta);
  else
    disp('running SARTA cloud MANY TIMES')
    prof = get_sarta_cloud100layerNtimes(h,ha,prof,pa,run_sarta);
  end
else
  disp('you did not ask for SARTA cloudy to be run; not changing p.rcalc')
end

tnow = toc;
fprintf(1,'TOTAL : %8.6f minutes to process \n',tnow/60);
