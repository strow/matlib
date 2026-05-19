function [h,prof] = do_subcol_cloudprofs(h,prof,xsubcol_frac);

%% do this ala PCRTM wrapper, so can loop over subcolumns

if h.ptype == 0
  good = find(prof.cc >= 0.001);
    prof.clwc(good) = prof.clwc(good) ./ prof.cc(good); 
    prof.ciwc(good) = prof.ciwc(good) ./ prof.cc(good); 

  bad = find(prof.cc < 0.001);
    prof.clwc(bad) = 0.0;
    prof.ciwc(bad) = 0.0;

  %% xsubcol_frac = 0 or 1
  prof.clwc = prof.clwc .* xsubcol_frac';
  prof.ciwc = prof.ciwc .* xsubcol_frac';

elseif h.ptype > 0
  good = find(prof.cc >= 0.001);
    prof.gas_201(good) = prof.gas_201(good) ./ prof.cc(good); 
    prof.gas_202(good) = prof.gas_202(good) ./ prof.cc(good); 

  bad = find(prof.cc < 0.001);
    prof.gas_201(bad) = 0.0;
    prof.gas_202(bad) = 0.0;

  [mmjunkA,nnjunkA] = size(prof.gas_201);
  [mmjunkX,nnjunkX] = size(xsubcol_frac);

  if mmjunkX == 101
    %% xsubcol_frac = 0 or 1
    prof.gas_201 = prof.gas_201 .* xsubcol_frac;
    prof.gas_202 = prof.gas_202 .* xsubcol_frac;
  else
    %% xsubcol_frac = 0 or 1
    prof.gas_201 = prof.gas_201 .* xsubcol_frac';
    prof.gas_202 = prof.gas_202 .* xsubcol_frac';
  end
  
end

