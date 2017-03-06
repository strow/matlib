function [prof,orig_slabs] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta)

%% modelled on MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE_CFRAC/readecmwf91_nearest_gasNcloud_slabprof.m
%% also see /asl/rtp_prod/airs/rtp/create_rcalc_ecm_cld_allfov.m
%
% function [prof,orig_slabs] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta)
%   takes an input [h ha p pa] which incudes cloud structure from (ERA/ECMWF) and
%   then runs the SARTA code
% input
%   h,ha,p,pa  = structures from rtp file
%   run_sarta  = optional arguments
% output
%   prof       = final profile struture, with cloud and clear rads as needed, and if needed, new cloud slab info
%              = this depends on ForceNewSlabs and whether there were cloud slab fields in inout "p" for SARTA cloudy
%  orig_slabs  = if there were cloud slab fields in input "p" structure, these are saved if there are new slab fields to be computed
%              = that way the user can then check things out and see if he/she wants to diff input slabs vs output slabs
%
% >>>>>>>>
% WARNING : before Sept 2015
%         :  when "p" originally came in from ERA/ECM, the geopphysical tcc field was set in p.cfrac, and
%            there were MANY p.cfrac of the order of 1
%         : so this "tcc = cfrac" info was used in "set_fracs_deffs.m" as seed for "fake_cfracs.m"
%           which means if cfrac is OVERWRITTEN before this you are sunk, as wrong info will go into fake_cfracs
%         : which is what happened : after running through driver_sarta_cloud_rtp for the first time, cloud slabs are done, 
%            >>> cfrac redone in fake_cfracs.m <<< and very few close to 1, and OUTPUT as CFRAC into the resulting output structure
%         : so if you input this re-jigged cloud "cfrac" into driver_sarta_cloud_rtp and asked for
%           another iteration of forcing NewCloudSlabs there was a notable difference in the seeds,
%           and so the window biases will change quite a bit (from 2 K to 8 K!!)
%         : IN OTHER WORDS if you wanted to test re-randomizing, it was always better to send in
%            pristine "p.cfrac" from ERA rather than re-jigged p.cfrac !!!!
% also if there were NO p,cfrac coming in originally, then get_orig_slabs_info tries to intiialize them
%   randomly, using size(p.stemp)  and end up with pretty bad biases!!!!
% in other words, incoming tcc from ERA/ECM (via cfrac in the input profile structure)
%   has quite a bit of useful info!!!!!!!!! that is used by
%          put_into_prof.m   reset_cloud_slab_with_cloud_profile.m  set_fracs_deffs.m  
% >>>>>>>>
% WARNING : after Sept 2015
%         : fill_era and fill_ecm will now output p.tcc INSTEAD of p.cfrac, as well as p.cc, p.ciwc, p.clwc
%         : this "tcc" info is now correctly used in "set_fracs_deffs.m" as the seeds for "fake_cfracs.m"
%         : ****** p1.cfrac+p1.cfrac2-p1.cfrac12 = p1.tcc **************
%         : ****** p1.cfrac+p1.cfrac2-p1.cfrac12 = p1.tcc **************
%         : ****** p1.cfrac+p1.cfrac2-p1.cfrac12 = p1.tcc **************
% >>>>>>>>
% WARNING : all random seeds must be done OUTSIDE the driver code : suggest rng('shuffle','twister')
% >>>>>>>>
%
% run_sarta = optional structure argument that says
% >>> options for SARTA runs
%     *** run_sarta.tcc ***         = +1 if this field exists, then reset p.cfrac with this for making slabs
%                                         if redoing slabs, this is THE most important variable!!!!
%     run_sarta.clear               = +/-1 for yes/no, results into prof.clearcalc (DEFAULT -1)
%     run_sarta.cloud               = +/-1 for yes/no, results into prof.rcalc     (DEFAULT +1)
%     run_sarta.cumsum = -1           : ORIG DEFAULT go with "ecmwf2sarta" results (default before March 2012), mean of ciwc/clwc "GeorgeAumann pick"
%                        0 -- 1        : set cloud pressure based on cumulative sum of p.ciwc and p.clwc, 
%                        >  1--9998    : go for where cumsum(cloudOD) ~ N/100 (if that can be found)
%                        >= +9999      : NEW DEFAULT go for peak of wgt fcn of cloud ice, cloud liquid (ie go HIGH in atm, good for DCC) "Strow pick"
%                        <= -9999      :             go for peak of wgt fcn of cloud ice, cloud liquid (ie go HIGH in atm, good for DCC)
%                                      :             differs from +9999 ice clouds are limited be 400 mb >= icetop >= 0 ("Strow pick" has 1000 mb > icetop > 0)
%     run_sarta.cfrac < 0              : use random (DEFAULT)
%                     > 0 to < 1       : use fixed amount specified by user
%     run_sarta.klayers_code        = string to klayers
%     run_sarta.sartaclear_code     = string to sarta clear executable
%     run_sarta.sartacloud_code     = string to sarta cloud executable
%     run_sarta.ice_water_separator = set all ciwc/clwc to ice above this, water below this 
%                                     (DEFAULT = -1, use ciwc/clwc structures as is)
%     run_sarta.randomCpsize        = +1 (DEFAULT) to randomize BOTH ice (based on Tcld) and water deff
%                                   = 20,   then water is ALWAYS 20 um (as in PCRTM wrapper),random ice
%                                           (based on Tcld)
%                                   = -1,   then water is MODIS dme, random ice (based on Tcld)
%                                   = 9999, then water is MODIS dme, random ice (based on KNLiou Tcld)
%     run_sarta.co2ppm              = -1 to use 370 + (yy-2002)*2.2) in pcrtm/sarta
%                                   = +x to use user input value     in pcrtm/sarta 
%                                   =  0 to use DEFAULT klayers = 385 (set in executable by Scott; equivalent to run_sarta.co2ppm = +385)
%                   this is done to keep it consistent with PCRTM
%                   however, also have to make sure this is only enabled if h.glist does NOT include gas_2
%      ForceNewSlabs                = -1 (default) to keep any slab clouds that are input, as they are
%                                   = +1           to force new slabs to be derived from clwc,ciwc,cc
%      Slab_or_100layer             = +1 for slab clouds/-1 for 100 layer clouds (which then need their own sarta,klayers,ncol,overlap)
%                                   
% >>> test ONE cloud
%     run_sarta.waterORice          = 0 (default, use both clouds)
%                                   = +1 turn off ice   clouds, keep water only, set ncol0 = -1, set cc = 1
%                                   = -1 turn off water clouds, keep ice   only, set ncol0 = -1, set cc = 1
%       if run_sarta.waterORice = +/-1 in PCRTM we set run_sarta.ncol0 == -1, p.cc = 1 and turn off water or ice clouds
%          while in SARTA it turns off appropriate ice or water slab
%          this is test of ONE SLAB CLOUD vs ONE COLUMN CLOUD
%
% Requirements : 
%   p must contain ciwc clwc cc from ERA/ECMWF (ie 91xN or 37xN) as well as gas_1 gas_3 ptemp etc
%   STRONGLY RECOMMENDED : p.cfrac from ERA/ECM is MUCH BETTER as a seed than what comes out from this routine
%   this code puts in its own particles sizes, cloud fracs and cloud tops based on ciwc,clwc,cc
%   h.ptype = 0 (ie must be levels profile)
%
% Can do arbitrary levels eg ECM (91 levels) or ERA (37 levels)
%
% can do "driver_sarta_cloud_rtp_onecldtest", but remember
%    simplest way of turing off ice   is set p.ciwc = 0
%    simplest way of turing off water is set p.clwc = 0,
% and then set p.cc = 1
%
% Written by Sergio DeSouza-Machado (with a lot of random cloud frac and dme by Scott Hannon)
%
% updates
%  03/06/2017 : introduce run_sarta.iNew_or_Orig_CXWC2OD : -1 (default) is to do OD = blah * qBlah / cc * diffZ                   almost PCRTM way
%                                                        : 0  is to do OD = blah * qBlah / cc * diffZ and then OD(cc < 1e-3) = 0  PCRTM way
%                                                        : +1 is to do OD = blah * qBlah * cc * diffZ                             SERGIO way
%  02/26/2017 : run_sarta.cumsum = N <= -9999 option allows the clouds to be placed at peak of cloud wgt fcn, with ice cloud top 400 mb >= icetop >= 0 mb
%  09/01/2015 : cleaned up random seeds (all must be done OUTSIDE the code), and include "tcc" as a run_sarta option
%  08/15/2015 : cleaned up code by including lots of common functions for driver_sarta_cloud_rtp.m driver_sarta_cloud100layer_rtp.m
%  08/14/2015 : f cprtop,cprbot,cngwat,cpsize exist, code no longer redoes profile--> slab unless forced to do so (run_sarta.ForceNewSlabs)
%  08/12/2015 : can send in CO2 and CH4 profiles (mostly from PCRTM code)
%  08/02/2015 : introduced "choose_klayers_sarta" to have default klayers/sarta
%
%  2014 : major fixes like fix problems with cfracs, cngwat
%  2014 : minor fixes like cut down on chatter, fix inf nan outputs
%
%  08/18/2013 : introduced run_sarta.randomCpsize, default = +1  to keep randomizing deff; 
%                                                             20 to keep ice = ice(T,pcrtm), water = 20 um
%                                                             -1 to use MODIS water DME
%  08/17/2013 : if cfrac does not exist, fix random initialization so as to make sure cfrac > 0 always
%  08/17/2013 : making more extensive use of run_sarta.cfrac by introducing it in set_fracs_deffs.m and check_for_errors.m
%  08/15/2013 : introduced run_sarta.ice_water_separator, so that everything above certain pressure == ice; 
%               below certain pressure = water
%  04/22/2013 : introduced run_sarta.cfrac so that we can eg have cfrac = 1
%  04/04/2013 : introduced/tweaked do_the_reset_cprtop_cloudOD.m
%  04/04/2013 : mega updates in check_for_errors.m
%  04/02/2013 : run_sarta.cumsum = N >= 9999 option allows the clouds to be placed at peak of cloud wgt fcn
%  03/28/2013 : run_sarta.cumsum = N > 1 option allows the clouds to be placed where cumsum(OD) ~ N/100 (so if N = 100, 
%               look for where cumsum(ODcld) ~ 1)
%  03/11/2013 : run_sarta.cumsum = 0--1 option allows the clouds to be placed where cumsum(ciwc)/sum(ciwc) = X
%                                  -1  option is to default wherever ecmwf2sarta put the cloud
%

%{
% testing the code
  [h,ha,p,pa] = rtpread('/asl/data/rtprod_airs/2012/05/01/cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp');
  [h,ha,p,pa] = rtpgrow(h,ha,p,pa);
  run_sarta.clear = +1;
  run_sarta.cloud = +1;
  run_sarta.cumsum = 100;   %% look for cumulativeOD = 100/100 = 1
  tic
  p1 = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
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
check_sarta_cloud_rtp_defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iAlreadyExistSlabClds = -1;  % assume no slab clouds
if isfield(p,'cprtop') & isfield(p,'cprbot') & isfield(p,'cpsize') & isfield(p,'cngwat') & isfield(p,'cfrac') & isfield(p,'ctype')
  disp('hmm : cprtop/cprbot,cpsize,cngwat,cfrac,ctype all exist')
  testA1 = p.cprtop; testA1(testA1 < 0) = NaN;
  testA2 = p.cprbot; testA2(testA2 < 0) = NaN;
  testA3 = p.cpsize; testA3(testA3 < 0) = NaN;
  testA4 = p.cngwat; testA4(testA4 < 0) = NaN;
  testA5 = p.cfrac;  testA3(testA3 < 0) = NaN;
  testA6 = p.ctype;  testA4(testA4 < 0) = NaN;

  testB1 = p.cprtop2; testB1(testB1 < 0) = NaN;
  testB2 = p.cprbot2; testB2(testB2 < 0) = NaN;
  testB3 = p.cpsize2; testB3(testB3 < 0) = NaN;
  testB4 = p.cngwat2; testB4(testB4 < 0) = NaN;
  testB5 = p.cfrac2;  testB3(testB3 < 0) = NaN;
  testB6 = p.ctype2;  testB4(testB4 < 0) = NaN;

  if nanmean(testA1) > 0 & nanmean(testA2) > 0 & nanmean(testA3) > 0 & nanmean(testA4) > 0 & nanmean(testA5) > 0 & nanmean(testA5) > 0
    disp('    and they are all valid non-negative means!');
    iAlreadyExistSlabClds = +1;
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pINPUT = p;

if run_sarta.ForceNewSlabs > 0
  %% force the making of new slab clouds
  if iAlreadyExistSlabClds > 0
    disp(' >>>> even though slab cloud params already exist, run_sarta.ForceNewSlabs > 0 ==> make new slab clouds')
  elseif iAlreadyExistSlabClds < 0
    disp(' >>>> slab cloud params do not exist, in any case run_sarta.ForceNewSlabs > 0 ==> make new slab clouds')
  end
  iAlreadyExistSlabClds = -1;
end

if iAlreadyExistSlabClds < 0
  %% need to add in slab cloud fields
  disp(' >>>>>>>>>>>>> adding in slab cloud fields <<<<<<<<<<<<<<<<<')
  [orig_slabs,p] = get_orig_slabs_info(p,run_sarta);
  prof = main_sarta_cloud_rtp(h,ha,p,pa,run_sarta,narginx); 
elseif iAlreadyExistSlabClds > 0
  %% slab cloud fields already exist, just run klayers and sarta
  disp(' >>>>>>>>>>>>> slab cloud fields already exist; simply running klayers and sarta <<<<<<<<<<<<<<<<<')
  orig_slabs = [];
  prof = p;
  main_compute_sarta_rads
end

tnow = toc;
fprintf(1,'TOTAL : %8.6f minutes to process \n',tnow/60);
