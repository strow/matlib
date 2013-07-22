function [head,prof] = ...
   level2slab_clouds_CFRAC(head,prof,rCfracType,rCNGWATiceMult,rCPSIZEiceMult)

%%[h,p] = level2slab_clouds_CFRAC(head,prof,[rCfracType,rCNGWATiceMult,rCPSIZEIceMult])
%  last three variables are optional
%
% INPUT
% head,  prof  : come from ERA or ECMWF reader
% rCfracType = +x;  %% use constant |x| for cfrac1,cfrac2, cfrac12 
%                    %% where x is between 0 and 1
%            = -1;  %% use wgtd avg of cc when doing cfrac,cfrac2,cfrac12
%            = -2;  %% use SAME wgtd avg of cc for all cfrac,cfrac2,cfrac12
%            = 0;   %% default, use Scott's methods {OPTIONAL; DEFAULT = 0}
% rCNGWATiceMult = {OPTIONAL; DEFAULT = 1} : multiplier for rCNGWAYIce
% rCPSIZEiceMult = {OPTIONAL; DEFAULT = 1} : multiplier for rCPSIZEIce
%
% used by eg RTPMAKE/loop_make_eracloudrtp_sergio_CFRAC.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
  rCfracType     = 0;       %% use Scott's method
  rCNGWATiceMult = 1.0;     %% optimal should be 10
  rCPSIZEiceMult = 1.0;     %% optimal should be 0.5
end
if nargin == 3
  rCNGWATiceMult = 1.0;     %% optimal should be 10
  rCPSIZEiceMult = 1.0;     %% optimal should be 0.5
end
if nargin == 4
  rCPSIZEiceMult = 1.0;     %% optimal should be 0.5
end

if rCfracType > 0 & rCfracType <= 1
  %% user wants fixed cfrac between 0 and 1
  rX = abs(rCfracType);
  if rX >= 0.99
    rX = 0.99;
  end
  if rX <= 0.01
    rX = 0.01;
  end
end

% Min allowed cloud fraction
cmin = 0.001;

% Max allowed cngwat[1,2]
cngwat_max = 500;

disp('****************')
disp('adding on clouds')
disp('****************')

%% basically copied from mkcldfile_era_hannon2.m

[nlev nprof] = size(prof.ptemp);

profX = prof;
ecmwfcld2sartacld
% wait quite a bit of time

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set rtpV201 cloud2 fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%
prof.cngwat2 = prof.udef(11,:);
prof.cpsize2 = prof.udef(12,:);  % replaced later
prof.cprtop2 = prof.udef(13,:);
prof.cprbot2 = prof.udef(14,:);
prof.cfrac2  = prof.udef(15,:);  % replaced later
prof.cfrac12 = prof.udef(16,:);  % replaced later
prof.ctype2  = prof.udef(17,:);

% Replace cfrac info
if rCfracType == 0
  %%%%%%%%%%%%%% default : Scott

  [cfracw, cfraci] = total_cfrac(profX.plevs,profX.cc,profX.clwc,profX.ciwc);
  tcc = profX.cfrac;

  [prof.cfrac, prof.cfrac2, prof.cfrac12] = ...
     fake_cfracs(tcc, cfracw, cfraci,prof.ctype, prof.ctype2);
  %%%%%%%%%%%%%% default : Scott

elseif rCfracType == -1
  % Replace cfrac info using CC
  [cfracw, cfraci] = total_cfrac(profX.plevs,profX.cc,profX.clwc,profX.ciwc);
  tcc = quick_cloud_frac(profX);

  [prof.cfrac, prof.cfrac2, prof.cfrac12] = ...
       fake_cfracs(tcc, cfracw, cfraci,prof.ctype, prof.ctype2);

elseif rCfracType == -2
  tcc = quick_cloud_frac(profX);
  prof.cfrac = tcc;
  prof.cfrac2 = tcc;
  prof.cfrac12 = tcc;

else
  %% assume constant value
  tcc = quick_cloud_frac(profX);
  prof.cfrac(tcc > 0)   = rX;
  prof.cfrac2(tcc > 0)  = rX;
  prof.cfrac12(tcc > 0) = rX;

  prof.cfrac(tcc <= 0)   = 0.0;
  prof.cfrac2(tcc <= 0)  = 0.0;
  prof.cfrac12(tcc <= 0) = 0.0;
end


clear profX

% Compute cloud temperature
pavg = 0.5*(prof.cprtop + prof.cprbot);
ibad = find(prof.cfrac == 0);
pavg(ibad) = 500; % safe dummy value
tavg1 = rtp_t_at_p(pavg, head, prof);

pavg = 0.5*(prof.cprtop2 + prof.cprbot2);
ibad = find(prof.cfrac2 == 0);
pavg(ibad) = 500; % safe dummy value
tavg2 = rtp_t_at_p(pavg, head, prof);
clear ibad


fprintf(1,'---->> rCNGWATiceMult = %8.6f rCPSIZEiceMult = %8.6f <<----\n',rCNGWATiceMult,rCPSIZEiceMult)

% Replace cpsize and cpsize2
iceflag = zeros(1,nprof);  %% zero === WATER
ii = find(prof.ctype == 201);
  iceflag(ii) = 1;           %% replace with 1 == ICE where needed
prof.cpsize  = fake_cpsize(tavg1, iceflag, 1);
  prof.cngwat(ii) = prof.cngwat(ii) * rCNGWATiceMult;
  prof.cpsize(ii) = prof.cpsize(ii) * rCPSIZEiceMult;

iceflag = zeros(1,nprof);  %% zero === WATER
ii = find(prof.ctype2 == 201);
  iceflag(ii) = 1;           %% replace with 1 == ICE where needed
prof.cpsize2 = fake_cpsize(tavg2, iceflag, 1);
  prof.cngwat2(ii) = prof.cngwat2(ii) * rCNGWATiceMult;
  prof.cpsize2(ii) = prof.cpsize2(ii) * rCPSIZEiceMult;

clear iceflag tavg1 tavg2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now do error checking

% Remove cloud fractions less than some minimum
hcmin = 0.5*cmin;

ix=find(prof.cfrac < cmin);
prof.cfrac(ix)  = 0;
prof.cfrac12(ix)= 0;
prof.cngwat(ix) = 0;
prof.ctype(ix)  = -1;
prof.cprtop(ix) = -9999;
prof.cprbot(ix) = -9999;

ix=find(prof.cfrac2 < cmin);
prof.cfrac2(ix)  = 0;
prof.cfrac12(ix) = 0;
prof.cngwat2(ix) = 0;
prof.ctype2(ix)  = -1;
prof.cprtop2(ix) = -9999;
prof.cprbot2(ix) = -9999;

ix = find(prof.cfrac12 >= hcmin & prof.cfrac12 < cmin);
prof.cfrac12(ix) = cmin;
ix = find(prof.cfrac12 < hcmin);
prof.cfrac12(ix) = 0;
junk = prof.cfrac(ix) + prof.cfrac2(ix);
ii = ix( find(junk > 1) );
ii1 = ii( find(prof.cfrac(ii) > prof.cfrac2(ii)) );
ii2 = setdiff(ii,ii1);
prof.cfrac(ii1) = prof.cfrac(ii1)-hcmin;
prof.cfrac2(ii2) = prof.cfrac2(ii2)-hcmin;

for ixx = 1 : 2
  % Error : CPRBO1 > CPRTO2
  prof.cprbot(prof.cprbot >= prof.cprtop2) = ...
    prof.cprtop2(prof.cprbot >= prof.cprtop2)- 50;

  % Error: CPRBOT > SPRES
  prof.cprbot(prof.cprbot >= prof.spres) = ...
    prof.spres(prof.cprbot  >= prof.spres)- 50;
  prof.cprbot2(prof.cprbot2 >= prof.spres) = ...
    prof.spres(prof.cprbot2 >= prof.spres)- 50;

  % Error: CPRTOP > CPRBOT
  prof.cprtop(prof.cprtop >= prof.cprbot) = ...
    prof.cprbot(prof.cprtop   >= prof.cprbot) - 50;
  prof.cprtop2(prof.cprtop2 >= prof.cprbot2) = ...
    prof.cprbot2(prof.cprtop2 >= prof.cprbot2)- 50;

  % now should not need these two
  % CPRTO1 outside allowed PLEVS1 to  SPRES range
  % prof.cprtop(prof.cprtop > prof.spres) = ...
  %   prof.spres(prof.cprtop > prof.spres);
  % extra check
  % prof.cprtop(prof.cprtop > prof.cprbot) = ...
  %   prof.cprbot(prof.cprtop > prof.cprbot);

  %% make sure everything at -9999
  prof.cprtop(prof.cprtop < -9999) = -9999;
  prof.cprbot(prof.cprbot < -9999) = -9999;
  prof.cprtop2(prof.cprtop2 < -9999) = -9999;
  prof.cprbot2(prof.cprbot2 < -9999) = -9999;
  prof.cfrac(prof.cprtop <= -9999) = 0.0;
  prof.cfrac(prof.cprbot <= -9999) = 0.0;
  prof.cfrac2(prof.cprtop2 <= -9999) = 0.0;
  prof.cfrac2(prof.cprbot2 <= -9999) = 0.0;

  %% make sure everything at -9999
  prof.cprtop(prof.cprtop < 0) = -9999;
  prof.cprbot(prof.cprbot < -9999) = -9999;
  prof.cprtop2(prof.cprtop2 < -9999) = -9999;
  prof.cprbot2(prof.cprbot2 < -9999) = -9999;
  prof.cfrac(prof.cprtop <= -9999) = 0.0;
  prof.cfrac(prof.cprbot <= -9999) = 0.0;
  prof.cfrac2(prof.cprtop2 <= -9999) = 0.0;
  prof.cfrac2(prof.cprbot2 <= -9999) = 0.0;
end
