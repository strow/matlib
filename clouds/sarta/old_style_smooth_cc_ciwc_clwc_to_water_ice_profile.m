function [watercld,icecld,aa,ptemp,cut440] = old_style_smooth_cc_ciwc_clwc_to_water_ice_profile(xcumsum,profX,ii,nlev,iSmooth,rGaussianCutoff);

%%% this has beeen superseded by new_style_smooth_cc_ciwc_clwc_to_water_ice_profile *********

iSmoothX = iSmooth;
blahW = profX.clwc(:,ii)/(max(profX.clwc(:,ii))+1e-16);
blahI = profX.ciwc(:,ii)/(max(profX.ciwc(:,ii))+1e-16);
numI = length(find(blahI > rGaussianCutoff*blahI));
numW = length(find(blahW > rGaussianCutoff*blahW));

iDo = -1;
if nlev > 20
  if numI > 3 & numW > 3
    iSmoothX = iSmooth;
    iDo = +1;
  elseif numW > 3 & numI == 0
    iSmoothX = iSmooth;
    iDo = +1;
  elseif numI > 3 & numW == 0
    iSmoothX = iSmooth;
    iDo = +1;
  else
    iSmoothX = 0.5;
    iSmoothX = 2;
    iDo = -1;         %% no ice or water
  end
end

if iDo > 0
   % lots of points
  [shiftedx,shiftedy,plevs]    = smooth1aa(1:nlev,profX.plevs(:,ii)',iSmoothX);
  [shiftedx,shiftedy,watercld] = smooth1aa(1:nlev,profX.clwc(:,ii)' ,iSmoothX);
  [shiftedx,shiftedy,icecld]   = smooth1aa(1:nlev,profX.ciwc(:,ii)' ,iSmoothX);

  %npoly = 1;    nframe = 3;
  %[sgplevs]    = sgolayfilt(double(profX.plevs(:,ii)') ,npoly,nframe);
  %[sgwatercld] = sgolayfilt(double(profX.clwc(:,ii)')  ,npoly,nframe);
  %[sgicecld]   = sgolayfilt(double(profX.ciwc(:,ii)')  ,npoly,nframe);
else
  % few points
  plevs    = profX.plevs(:,ii);
  watercld = profX.clwc(:,ii);
  icecld   = profX.ciwc(:,ii);
end

aa        = [];
ptemp     = [];
cutoff440 = [];