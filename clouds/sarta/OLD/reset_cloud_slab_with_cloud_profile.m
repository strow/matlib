function [h,prof] = reset_cloud_slab_with_cloud_profile(h,prof,xfrac);

%% this function takes the info in the SLAB profile, and replaces it with 
%% a PROFILE based on ciwc, clwc ... very simple
%%
%% only smart thing about this is that if xfrac >= 0, it replaces the cloud fraction with
%% xfrac. That way you can eg set xfrac = 1, and then do subcolumn simulations

cloud_gunit = [21 21]';
h.ngas  = h.ngas + 2;
h.glist = [h.glist;  [201 202]'];
h.gunit = [h.gunit;  cloud_gunit];

prof.ctype   = ones(size(prof.stemp)) * 201;
prof.ctype2  = ones(size(prof.stemp)) * 101;
prof.gas_201 = prof.ciwc;  sum201 = nansum(prof.ciwc);
prof.gas_202 = prof.clwc;  sum202 = nansum(prof.clwc);

if isfield(prof,'cprtop')
  prof = rmfield(prof,'cprtop');
end

if isfield(prof,'cprtop2')
  prof = rmfield(prof,'cprtop2');
end

if isfield(prof,'cprbot')
  prof = rmfield(prof,'cprbot');
end

if isfield(prof,'cprbot2')
  prof = rmfield(prof,'cprbot2');
end

if xfrac >= 0
  prof.cfrac  = zeros(size(prof.cfrac));  oo = find(sum201 > 0); prof.cfrac(oo) = xfrac;
  prof.cfrac2 = zeros(size(prof.cfrac2)); oo = find(sum202 > 0); prof.cfrac2(oo) = xfrac;
  prof.cfrac12 = max(prof.cfrac,prof.cfrac2);
  oo = find(sum201 <= 0); prof.cfrac12(oo) = prof.cfrac(oo);
  oo = find(sum202 <= 0); prof.cfrac12(oo) = prof.cfrac2(oo);

  %oo = find(prof.cfrac  > 0 & prof.cngwat > 0);  prof.cfrac(oo)  = xfrac;
  %oo = find(prof.cfrac2 > 0 & prof.cngwat2 > 0); prof.cfrac2(oo) = xfrac;
  %oo = find(prof.cfrac  > 0 & prof.cfrac2 > 0);  prof.cfrac12(oo)  = xfrac;
end
