function prof = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta)

%% modelled on MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE_CFRAC/readecmwf91_nearest_gasNcloud_slabprof.m
%% also see /asl/rtp_prod/airs/rtp/create_rcalc_ecm_cld_allfov.m

% function prof = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta)
% takes an input [h ha p pa] which incudes cloud structure from (ERA/ECMWF) and
% then runs the SARTA code
%
% run_sarta = optional structure argument that says
%   run_sarta.clear = +/-1 for yes/no
%   run_sarta.cloud = +/-1 for yes/no
%   run_sarta.klayers_code = string to klayers
%   run_sarta.sartaclear_code = string to sarta clear executable
%   run_sarta.sartacloud_code = string to sarta cloud executable
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
  [h,ha,p,pa] = oldrtpread('/asl/data/rtprod_airs/2012/05/01/cld_ecm_41ch.airs_ctr.2012.05.01.10.rtp');
  [h,ha,p,pa] = rtpgrow(h,ha,p,pa);
  run_sarta.cloud = +1;
  tic
  p1 = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
  toc
  rtpwrite('/asl/data/rtprod_airs/2012/05/01/pcrtm_cld_ecm_41ch.airs_ctr.2012.05.01.10_sarta.rtp',h,ha,p1,pa);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defaults
if nargin == 4
  run_sarta.clear = -1;
  run_sarta.cloud = -1;
elseif nargin == 5
  if ~isfield(run_sarta,'clear')
    run_sarta.clear = -1;
  end
  if ~isfield(run_sarta,'cloud')
    run_sarta.cloud = -1;
   end
  if ~isfield(run_sarta,'klayers_code')
    run_sarta.klayers_code = '/asl/packages/klayers/Bin/klayers_airs'; 
  end   
  if ~isfield(run_sarta,'sartaclear_code')
    run_sarta.sartaclear_code = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';
  end
  if ~isfield(run_sarta,'sartacloud_code')
    run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
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
prof  = p;  %% this is the output  <<<<<<<<<<<<<<<<<<<<<

profX = p;  %% temporary

nlev     = ceil(mean(profX.nlevs));
nlev_std = (std(double(profX.nlevs)));

if nlev_std > 1e-3
  error('oops : code assumes ERA (37 levs) or ECMWF (91 levs) or other constant numlevs model')
end

if h.ptype ~= 0
  error('need levels input!')
end

tic
ecmwfcld2sartacld;        %% figure the two slab cloud profile info here, using profX
                          %% this then puts the info into "prof" by calling put_into_prof
  prof.ctype  = double(prof.ctype);
  prof.ctype2 = double(prof.ctype2);

put_into_V201cld_fields   %% puts cloud info from above into rtpv201 fields 
set_fracs_deffs           %% sets fracs and particle effective sizes eg cfrac2, and clears profX
check_for_errors          %% see if there are possible pitfalls

if run_sarta.cloud > 0
  disp('running SARTA cloud')
  get_sarta_cloud;
  prof.rcalc = profRX2.rcalc;
else
  disp('you did not ask for SARTA cloudy to be run; not changing p.rcalc')
end

tnow = toc;
fprintf(1,'TOTAL : %8.6f minutes to process \n',tnow/60);

 
