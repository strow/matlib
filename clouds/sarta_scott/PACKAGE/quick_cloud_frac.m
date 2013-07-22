function [ccx,cfracx] = quick_cloud_frac(p0);

%% TCC is total cloud cover (1 x nprofs),     gets stuffed into p.cfrac
%%  CC is cloud cover,       (nlevs x nprofs) field, used in
%%                           level2slab_clouds_CFRAC

%% ccx uses cc :     ccx = 1 - prod(i=1:nlevs) (1-cc(i))
%%   so eg if cloud fraction cc(i) = 0 at all layers, ccx = 1-prod(1.1.1) = 0
%%   so eg if cloud fraction cc(i) = 1 at all layers, ccx = 1-prof(0.0.0) = 1
cfracx = p0.cfrac;

[mm,nn] = size(p0.cc);
ccx = p0.cc;
ccx = 1-ccx;
ccx = 1 - prod(ccx);

plot(1:nn,ccx,1:nn,cfracx)