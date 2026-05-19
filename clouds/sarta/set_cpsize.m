function prof = set_cpsize(head, profx,randomCpsize)

%% randomCpsize = 1    ==> random dme_water (centered about 20 um), dme_ice = based on Scott Tcld, randomized
%%                20   ==> dme_water FIXED at 20 um, dme_ice = based on KNLiou Tcld, randomized
%%                -1   ==> dme_water based on MODIS climatology, dme_ice = based on Scott Tcld, randomized
%%              9999   ==> dme water based on MODIS climatology, dme_ice = based on KNLiou Tcld
%% note this in now INDPT of runx.cumsum ie you placed the cprtop/cprbot according to runx.cumsum ... now code is finding the deff

prof = profx;

% Compute cloud temperature
pavg = 0.5*(prof.cprtop + prof.cprbot);
ibad = find(prof.cfrac == 0);
pavg(ibad) = 500; % safe dummy value
tavg1 = rtp_t_at_p_sarta2slab(pavg, head, prof);

pavg = 0.5*(prof.cprtop2 + prof.cprbot2);
ibad = find(prof.cfrac2 == 0);
pavg(ibad) = 500; % safe dummy value
tavg2 = rtp_t_at_p_sarta2slab(pavg, head, prof);

clear ibad

%%%% >>>>>>>>>>>>>>>>>>>>>>>>> guts of cpsize setting >>>>>>>>>>>>>>>>>>>>>>>>>
nprof = length(prof.stemp);

% Replace cpsize and cpsize2
iceflag = zeros(1,nprof);
ii = find(prof.ctype == 201);
iceflag(ii) = 1;
prof.cpsize  = fake_cpsize(tavg1, iceflag, randomCpsize);

iceflag = zeros(1,nprof);
ii = find(prof.ctype2 == 201);
iceflag(ii) = 1;
prof.cpsize2 = fake_cpsize(tavg2, iceflag, randomCpsize);
clear iceflag tavg1 tavg2

if (randomCpsize == -1) | (randomCpsize == +9999)
  disp('>>>>>>> resetting water dme according to water climatology')
  prof = modisL3_map_rtp(prof);
end
%%%% >>>>>>>>>>>>>>>>>>>>>>>>> guts of cpsize setting >>>>>>>>>>>>>>>>>>>>>>>>>
