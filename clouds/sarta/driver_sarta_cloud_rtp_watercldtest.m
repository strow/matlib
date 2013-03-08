function [prof,index] = driver_sarta_cloud_rtp_watercldtest(h,ha,p,pa,run_sarta,waterORice)

%% modelled on MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE_CFRAC/readecmwf91_nearest_gasNcloud_slabprof.m
%% also see /asl/rtp_prod/airs/rtp/create_rcalc_ecm_cld_allfov.m

% function prof = driver_sarta_cloud_rtp_watercldtest(h,ha,p,pa,run_sarta,waterORice)
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
%
% Requirements : 
%   p must contain ciwc clwc cc from ERA/ECMWF (ie 91xN or 37xN) as well as gas_1 gas_3 ptemp etc
%   this code puts in its own particles sizes, cloud fracs and cloud tops based on ciwc,clwc,cc
%   h.ptype = 0 (ie must be levels profile)
%
% Can do ECM (91 levels) or ERA (37 levels)
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

  iAIRS = +1;
  iAIRS = -1;
 
  if iAIRS > 0
    fname = '/asl/data/rtprod_airs/2012/05/01/cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp';
    [h,ha,p,pa] = rtpread(fname);
    [h,ha,p,pa] = rtpgrow(h,ha,p,pa);
  else
    fname = '/strowdata1/shared/imbiriba/test_file_for_sergio.rtp';
    run_sarta.sartaclear_code = '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_wcon_nte';
    run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
    [h,ha,p,pa] = rtpread(fname);
  end

  [h,p] = subset_rtp_allcloudfields(h,p,[],[],1:415);

  run_sarta.cloud = +1;
  run_sarta.clear = +1;
  if iAIRS > 0
    [px,ix] = driver_sarta_cloud_rtp_watercldtest(h,ha,p,pa,run_sarta);

    run_sarta.clear = +1;
    run_sarta.cloud = +1;
    addpath ../pcrtm
    addpath ../sarta
    cd ../pcrtm
    p1 = driver_pcrtm_cloud_rtp(h,ha,px,pa,run_sarta);
    cd ../sarta

    p1.rsarta_water_cloud_only = px.rcalc;

    figure(1)
    plot(h.vchan,rad2bt(h.vchan,p1.rclearcalc)-rad2bt(h.vchan,p1.rcalc))
    figure(2)
    plot(h.vchan,rad2bt(h.vchan,p1.rclearcalc)-rad2bt(h.vchan,p1.rsarta_water_cloud_only))
    title('cloud effect')
    figure(3)
    %% rcalc is input, and NOT changed by driver_pcrtm_cloud_rtp
    plot(h.vchan,rad2bt(h.vchan,p1.rcalc)-rad2bt(h.vchan,p1.rsarta_water_cloud_only))  
    %% better be the same, both done by SARTA
    plot(h.vchan,rad2bt(h.vchan,p1.rclearcalc)-rad2bt(h.vchan,p1.sarta_clear))
    %% better be close, as this is clearsky
    plot(h.vchan,rad2bt(h.vchan,p1.rad_clrsky)-rad2bt(h.vchan,p1.sarta_clear))

    %% better be zero, as this is cloudy done by SARTA
    plot(h.vchan,rad2bt(h.vchan,p1.sarta_cloud)-rad2bt(h.vchan,p1.rsarta_water_cloud_only))

    figure(4)
    plot(h.vchan,rad2bt(h.vchan,p1.rad_allsky)-rad2bt(h.vchan,p1.sarta_cloud))
    title('PCRTM watercld - SARTA watercld')

  else
    [px] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
  end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defaults
if nargin == 4
  run_sarta.clear = -1;
  run_sarta.cloud = +1;
  run_sarta.cumsum = -1;
  %run_sarta.klayers_code = '/asl/packages/klayers/Bin/klayers_airs';
  %run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  waterORice = +1; % keep only water clds
elseif nargin == 5
  waterORice = +1; % keep only water clds
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

prof = put_into_V201cld_fields(prof);    %% puts cloud info from above into rtpv201 fields 
  prof.ctype  = double(prof.ctype);
  prof.ctype2 = double(prof.ctype2);

prof = set_fracs_deffs(head,prof,profX,...
            cmin,cngwat_max);            %% sets fracs and particle effective sizes eg cfrac2

if run_sarta.cumsum > 0
  prof = reset_cprtop(prof);
end

prof = check_for_errors(prof);           %% see if there are possible pitfalls

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

 
