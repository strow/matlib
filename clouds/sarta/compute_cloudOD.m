function p1 = compute_cloudOD(p0,airslevels,airsheights,iNew_or_Orig_CXWC2OD);

if nargin == 3
  iNew_or_Orig_CXWC2OD =  0;  %%% change to OD = blah * qBlah / cc * diffZ; OD(cc < 1e-3) = 0 WHAT PCRTM DOES
  iNew_or_Orig_CXWC2OD = +1;  %%% change to OD = blah * qBlah * cc * diffZ                    Mar 2017 SERGIO
  iNew_or_Orig_CXWC2OD = -1;  %%% stick  to OD = blah * qBlah / cc * diffZ                    Pre March 2017  DEFAULT
end

disp('computing level ODs')
%% basically same as reset_cprtop_cloudOD.m
%% first part is same as reset_cprtop_clouddOD.m but then this routine does much more

%% computes cloud ODs based on formulas given by Xianglei and Xiuhong
%% see PCRTM_compute_for_AIRS_spectra.m
%%
%%
%%

p1 = p0;
p1.orig_ctop  = p1.cprtop;
p1.orig_ctop2 = p1.cprtop2;

for ii = 1 : length(p0.stemp)
  [p1,iceOD,waterOD] = ice_water_deff_od(p1,airslevels,airsheights,ii,iNew_or_Orig_CXWC2OD);
end

