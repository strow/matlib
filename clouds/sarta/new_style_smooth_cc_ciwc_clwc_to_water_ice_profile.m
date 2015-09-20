function [watercld,icecld,aa,ptemp,cut440] = new_style_smooth_cc_ciwc_clwc_to_water_ice_profile(xcumsum,profX,ii);

%%% old_style_smooth_cc_ciwc_clwc_to_water_ice_profile been superseded by this new code

%% slabs can be resolved better if there are more points
plevs = profX.plevs(:,ii);
if length(plevs < 80)
  plevsX = (plevs(1:end-1) + plevs(2:end))/2;
  plevs = sort([plevs; plevsX]);
end
if length(plevs < 80)
  plevsX = (plevs(1:end-1) + plevs(2:end))/2;
  plevs = sort([plevs; plevsX]);
end

if plevs(1) < plevs(10)
  cut440 = find(plevs >= 440,1);
else
  cut440 = find(plevs <= 440,1);
end

watercld = interp1(log10(profX.plevs(:,ii)),profX.clwc(:,ii),log10(plevs));
icecld   = interp1(log10(profX.plevs(:,ii)),profX.ciwc(:,ii),log10(plevs));
ptemp    = interp1(log10(profX.plevs(:,ii)),profX.ptemp(:,ii),log10(plevs));

if ~exist('aa','var')
  aa = [];
end
aa = cloud_mean_press(aa,xcumsum,icecld,watercld,plevs,ii);

