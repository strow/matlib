function prof = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta)

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
%     run_sarta.cumsum = < 0        : just go with "ecmwf2sarta" results (default before March 2012)
%                        0 -- 1     : set cloud pressure based on cumulative sum of p.ciwc and p.clwc, 
%                        >  1--9998 : go for where cumsum(cloudOD) ~ N/100 (if that can be found)
%                        >= 9999    : go for peak of wgt fcn of cloud ice, cloud liquid
%     run_sarta.klayers_code = string to klayers
%     run_sarta.sartaclear_code = string to sarta clear executable
%     run_sarta.sartacloud_code = string to sarta cloud executable
%
% Requirements : 
%   p must contain ciwc clwc cc from ERA/ECMWF (ie 91xN or 37xN) as well as gas_1 gas_3 ptemp etc
%   this code puts in its own particles sizes, cloud fracs and cloud tops based on ciwc,clwc,cc
%   h.ptype = 0 (ie must be levels profile)
%
% Can do arbitrary levels eg ECM (91 levels) or ERA (37 levels)
%
% Written by Sergio DeSouza-Machado (with a lot of random cloud frac and dme by Scott Hannon)
%
% updates
%  04/04/2013 : introduced/tweaked do_the_reset_cprtop_cloudOD.m
%  04/04/2013 : mega updates in check_for_errors.m
%  04/02/2013 : run_sarta.cumsum = N >= 9999 option allows the clouds to be placed at peak of cloud wgt fcn
%  03/28/2013 : run_sarta.cumsum = N > 1 option allows the clouds to be placed where cumsum(OD) ~ N/100 (so if N = 100, look for where cumsum(ODcld) ~ 1)
%  03/11/2013 : run_sarta.cumsum = 0--1 option allows the clouds to be placed where cumsum(ciwc)/sum(ciwc) = X
%                                  -1  option is to default wherever ecmwf2sarta put the cloud
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defaults
if nargin == 4
  %% default to running sarta_cloudy
  run_sarta.clear = -1;
  run_sarta.cloud = +1;
  run_sarta.cumsum = -1;
  %run_sarta.klayers_code = '/asl/packages/klayers/Bin/klayers_airs';
  %run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  run_sarta.sartacloud_code = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
elseif nargin == 5
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
cmin = 0.001;

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

tic
[prof,profX] = ecmwfcld2sartacld(p,nlev,run_sarta.cumsum);   %% figure the two slab cloud profile info here, using profX
                                            %% this then puts the info into "prof" by calling put_into_prof w/in routine

%disp('check 1')
%junky = 199;
%junky
%[prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky)]
%disp('check 1')

prof = put_into_V201cld_fields(prof);    %% puts cloud info from above into rtpv201 fields 
  prof.ctype  = double(prof.ctype);
  prof.ctype2 = double(prof.ctype2);

%disp('check 2')
%[prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky)]
%disp('check 2')

prof = set_fracs_deffs(head,prof,profX,...
            cmin,cngwat_max);            %% sets fracs and particle effective sizes eg cfrac2

%disp('check 3')
%[prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky)]
%disp('check 3')

if run_sarta.cumsum > 0 & run_sarta.cumsum <= 1
  %% set cloud top according to cumulative sum fraction of ciwc or clwc
  profXYZ = prof;
  prof = reset_cprtop(prof);
elseif run_sarta.cumsum > 1
  %% set cloud top according to where cumulative cloud OD = run_sarta.cumsum/100, 
  %%      or if >= 9999, set at peak of cloud wgt fcn
  profXYZ = prof;
  prof = reset_cprtop_cloudOD(prof,run_sarta.cumsum/100);  
%{

  disp(' ')
  [prof.cngwat(junky) prof.cngwat2(junky)]
  [prof.cfrac(junky) prof.cfrac2(junky) prof.cfrac12(junky)]
  [prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky)]
  [prof.sarta_wgtpeakI(junky) prof.sarta_wgtpeakW(junky)]
  [profXYZ.cprtop(junky) profXYZ.cprbot(junky) profXYZ.cprtop2(junky) profXYZ.cprbot2(junky)]
  keyboard

  aaaa = load('/home/sergio/test_paul_p2.mat');
  figure(1); plot(prof.sarta_lvlODice,'b');   hold on; plot(aaaa.p2.pcrtm_lvlODice,'r'); hold off
  figure(2); plot(prof.sarta_lvlODwater,'b'); hold on; plot(aaaa.p2.pcrtm_lvlODwater,'r'); hold off
  ice = find(prof.ctype == 201);
  water = find(prof.ctype2 == 101);
  figure(3); plot(prof.sarta_lvl_iceOD_1(ice),prof.cprtop(ice),'bo',...
                  prof.sarta_lvl_waterOD_1(water),prof.cprtop2(water),'rx',...
                  [0 1000],[00 1000],'k')
             axis([0 1000 0 1000]); grid
  xlabel('OD = 1'); ylabel('cprtop')
  keyboard
%}
end

disp('---> checking cprtop vs cprbot vs spres')
iNotOK = +1;
iFix = 0;
  %[prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky)]
while iFix < 10 & iNotOK > 0
  iFix = iFix + 1;
  [prof,iNotOK] = check_for_errors(prof);           %% see if there are possible pitfalls x1
  %[prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky)]
  fprintf(1,' did n=%2i try at checking clouds \n',iFix)
end
if iFix > 10 & iNotOK > 0
  error('oops, could not fix cprtop vs cprbot vs spres')
end

clear profX

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  tic
  get_sarta_clear;
  toc
  prof.sarta_rclearcalc = profRX2.rcalc;
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

 
