function [h,prof] = do_subcol_cloudprofs(h,prof,xsubcol_frac);

%% do this ala PCRTM wrapper, so can loop over subcolumns

good = find(prof.cc >= 0.001);
  prof.clwc(good) = prof.clwc(good) ./ prof.cc(good); 
  prof.ciwc(good) = prof.ciwc(good) ./ prof.cc(good); 

bad = find(prof.cc < 0.001);
  prof.clwc(bad) = 0.0;
  prof.ciwc(bad) = 0.0;

%% xsubcol_frac = 0 or 1
prof.clwc = prof.clwc .* xsubcol_frac';
prof.ciwc = prof.ciwc .* xsubcol_frac';
