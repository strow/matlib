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
%        (default = -1, use ciwc/clwc structures as is)
%     run_sarta.randomCpsize        = +1 or 0 to turn on/off randomizing of ice and water deff
%                                      if 0, then water is ALWAYS 20 um (as in PCRTM wrapper)
%
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
%  08/18/2013 : introduced run_sarta.randomCpsize, default = +1 to keep randomizing deff; 
%                                                             0 to keep ice = ice(T,pcrtm), water = 20 um
%  08/17/2013 : if cfrac does not exist, fix random initialization so as to make sure cfrac > 0 always
%  08/17/2013 : making more extensive use of run_sarta.cfrac by introducing it in set_fracs_deffs.m and check_for_errors.m
%  08/15/2013 : introduced run_sarta.ice_water_separator, so that everything above certain pressure == ice; 
%               below certain pressure = water
%  04/22/2013 : introduced run_sarta.cfrac so that we can eg have cfrac = 1
%  04/04/2013 : introduced/tweaked do_the_reset_cprtop_cloudOD.m
%  04/04/2013 : mega updates in check_for_errors.m
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
% testing the code for 100layer clouds, random overlap
  [h,ha,p,pa] = rtpread('/asl/data/rtprod_airs/2012/05/01/cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp');
  [h,ha,p,pa] = rtpgrow(h,ha,p,pa);
  run_sarta.clear = +1;
  run_sarta.cloud = +1;
  run_sarta.ncol  = 50;
  run_sarta.cfrac = 1;
  run_sarta.ice_water_separator = 440;
  tic
  p1 = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
  toc
  rtpwrite('/asl/data/rtprod_airs/2012/05/01/pcrtm_cld_ecm_41ch.airs_ctr.2012.05.01.10_sarta.rtp',h,ha,p1,pa);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defaults
run_sarta.cumsum              = -1;  %% use pre-2012 cloudtop heights, without adjustments
  %% need this just to do basic calls for the code

if nargin == 4
  %% default to running sarta_cloudy100layer for 1 subcolumns
  run_sarta.clear               = -1;  %% do not run clear code
  run_sarta.cloud               = +1;  %% run cloudy code
  run_sarta.ice_water_separator = -1;  %% do not separate out ciwc and clwc by pressure; ie believe the NWP are correct
  run_sarta.randomCpsize        = +1;  %% keep randomizing dme for ice and water
  run_sarta.cfrac               = -1;  %% use random (instead of fixed) cfracs
  run_sarta.ncol                =  1;
  run_sarta.overlap             = +3;  %% maximal random overlap

%% SLAB
  %run_sarta.klayers_code    = '/asl/packages/klayers/Bin/klayers_airs';
  %run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  %run_sarta.klayers_code    = '/asl/packages/klayersV205/BinV201/klayers_airs';
  %run_sarta.sartacloud_code = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
%% SLAB

%% PROFILE
  %% see /home/sergio/klayersV205/Src_rtpV201_100layercloudamountsize
  run_sarta.klayers_code = ...
    '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';

  %% see /home/sergio/SARTA_CLOUDY/v108_Src_rtpV201_pclsam_slabcloud_hg3_100layerNEW
  run_sarta.sartacloud_code = ...
    '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';

elseif nargin == 5
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
  %if ~isfield(run_sarta,'cumsum')
  %  run_sarta.cumsum = -1;
  %end
  if ~isfield(run_sarta,'cfrac')
    run_sarta.cfrac = -1;
  end
  if ~isfield(run_sarta,'ncol')
    run_sarta.ncol = +1;
  end
  if ~isfield(run_sarta,'overlap')
    run_sarta.overlap = +3;
  end
  if ~isfield(run_sarta,'klayers_code')
    run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
    run_sarta.klayers_code = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';
  end   
  if ~isfield(run_sarta,'sartaclear_code')
    run_sarta.sartaclear_code = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';
  end
  if ~isfield(run_sarta,'sartacloud_code')
    run_sarta.sartacloud_code = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
    run_sarta.sartacloud_code = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';
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
  disp('>>>>>>>> warning : setting SEPARATOR for ice and water .... initializing')
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

tic
%% run_sarta.cumsum == -1
[prof,profX] = ecmwfcld2sartacld(p,nlev,run_sarta.cumsum);   %% figure the two slab cloud 
                  %% profile info here, using profX
                  %% this then puts the info into "prof" by calling put_into_prof w/in routine

prof = put_into_V201cld_fields(prof);    %% puts cloud info from above into rtpv201 fields 
  prof.ctype  = double(prof.ctype);
  prof.ctype2 = double(prof.ctype2);

%% sets fracs and particle effective sizes eg cfrac2
prof = set_fracs_deffs(head,prof,profX,cmin,cngwat_max,run_sarta.cfrac,run_sarta.randomCpsize);

%disp('3')
%[prof.cngwat]

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

disp('---> checking cprtop vs cprbot vs spres')
iNotOK = +1;
iFix = 0;
while iFix < 10 & iNotOK > 0
  iFix = iFix + 1;
  [prof,iNotOK] = check_for_errors(prof,run_sarta.cfrac);           %% see if there are possible pitfalls x1
  %[prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky)]
  fprintf(1,' did n=%2i try at checking clouds \n',iFix)
end
if iFix > 10 & iNotOK > 0
  error('oops, could not fix cprtop vs cprbot vs spres')
end

clear profX

%%%%%%%%%%%%%%%%%%%%%%%%%
%% now essentially IGNORE the above and instead do the cloud PROFILES stuff
[h,prof] = reset_cloud_slab_with_cloud_profile(h,prof,run_sarta.cfrac);

%% and add on aux info, such as OD etc
prof = cloudOD_for100layer(prof,run_sarta.cumsum/100);

%%%%%%%%%%%%%%%%%%%%%%%%%

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  tic
  get_sarta_clear;
  toc
  prof.sarta_rclearcalc = profRX2.rcalc;
end

if run_sarta.cloud > 0 
  if run_sarta.ncol == 1
    %% all you have to do is run ONCE
    disp('running SARTA cloud ONCE')
    tic
      get_sarta_cloud100layer;
    toc
    prof.rcalc = profRX2.rcalc;
  else
    profNCOL_IP = prof;   %%%% <<<<<<<<<<<<<<<< save this, need it a lot!!!!!
    h_IP        = h;
    get_sarta_cloud100layer_klayersONLY
    profNCOL_OP = pjunk;
    h_OP        = hjunk;
    [mmjunk,nnjunk] = size(profNCOL_OP.plevs);
    [unique_col_frac,ucol_num,ucol_num_same,subcol_frac] = ...
       get_subcolumn_frac_v2(length(prof.stemp), mmjunk, run_sarta.ncol, profNCOL_OP.cc',...
                                    run_sarta.overlap);

    for iCol = 1 : run_sarta.ncol
      fprintf(1,'running SARTA cloud N times : subcol %3i out of %3i \n',iCol,run_sarta.ncol);
      prof = profNCOL_OP;
      [hX,profX] = do_subcol_cloudprofs(h_OP,prof,squeeze(subcol_frac(:,iCol,:)));
      if hX.ptype == 0
        xciwc(iCol,:,:) = profX.ciwc;
        xclwc(iCol,:,:) = profX.clwc;
      else
        xciwc(iCol,:,:) = profX.gas_201;
        xclwc(iCol,:,:) = profX.gas_202;
      end
      get_sarta_cloud100layer_sartaONLY;
      %% junkcalc(iCol,:,:) = profRX2.rcalc;    %% MEMORY HOG
      if iCol == 1
        %% slower, but more memory efficient
        [sumy,sumysqr,Nmatr] = accum_mean_std(0,0,0,profRX2.rcalc,1);
      else
        %% slower, but more memory efficient
        [sumy,sumysqr,Nmatr] = accum_mean_std(sumy,sumysqr,Nmatr,profRX2.rcalc,iCol);
      end
    end

    prof = profNCOL_IP;

    %% prof.rcalc     = squeeze(nanmean(junkcalc,1));
    %% prof.rcalc_std = squeeze(nanstd(junkcalc,1));

    prof.rcalc = sumy./Nmatr;
    junk_mean = prof.rcalc;
    %prof.rcalc_std = ...
    %  real(sqrt((sumysqr - 2*junk_mean.*sumy + Nmatr.*junk_mean.*junk_mean)./(Nmatr-1)));
    prof.rcalc_std  = ...
      real(sqrt((sumysqr - 2*junk_mean.*sumy + Nmatr.*junk_mean.*junk_mean)./(Nmatr-0)));

  end
else
  disp('you did not ask for SARTA cloudy to be run; not changing p.rcalc')
end

tnow = toc;
fprintf(1,'TOTAL : %8.6f minutes to process \n',tnow/60);

 
