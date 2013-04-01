function [prof,profX] = ecmwfcld2sartacld(profIN,nlev,xcumsum);

%% called by readecmwf91_grid/nearest_gasNcloud.m
%%     "nlev" is set by readecmwf91_grid/nearest_gasNcloud
%%
%%%%%%%%%%%%%% will fail if iN or iW >= 6 at the smoothing  %%%%%%%%%%%%%%%%
%% first smooth plevs (P) , ice cloud profile (I) , water cloud profile (W)
%%   using iSmooth = 2
%% if num(maxN) = 5 or 6 after smoothing, boxshape dumps lowest M maxima
%%   so that there are 4 maxima at most
%% then take smoothed profiles and make them into boxshapes (N <= 4)
%% smooth 4 to 3 : if iN > 3 combines 4 clouds to 3
%%               : if iW > 3 combines 4 clouds to 3
%% smooth 3 to 2 : if iN > 2 combines 3 clouds to 2
%%               : if iW > 2 combines 3 clouds to 2
%% if we have total iN + iW > 2 then smooth 2 ice   clouds to 1
%%                                   smooth 2 water clouds to 1
%% finally combine the ice and water clouds
%%%%%%%%%%%%%% will fail if iN or iW >= 6 at the smoothing  %%%%%%%%%%%%%%%%
%% once it has changed the cloud profile into slabs, it calls "put_into_prof"
%% in that routine, the effective particle radii are preset .... eventually 
%% may make ice particle sizes vary with cloud top
%% ice_dme   = 60;
%% water_dme = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

profX = profIN;  %% dummy copy
prof  = profIN;  %% prof gets updated profile by profile, using "put_into_prof"

iLevsVers = 1;   %% orig, also has slight mistake in "put_into_profs"
iLevsVers = 2;   %% Oct 2011

if iLevsVers == 1
  if nlev > 20
    iSmooth = 2;
  else
    iSmooth = 1;
  end
elseif iLevsVers == 2
  iSmooth = 1;
end

rGaussianCutoff = 0.250;   %% done for ages, till July 2012; might put cloud too high
rGaussianCutoff = 0.375;   %% done after July 2012, trying to put cloud little lower
                           %% and therefore closer to maxpart of cloud

global iDoPlot
iDoPlot = +1;   %% plot stuff .. ugh slow
iDoPlot = -1;   %% do not plot stuff

iPrint = +1;    %% do     print chirpy talky comments
iPrint = -1;    %% do not print chirpy talky comments
tic;

jj = 0;
iiii = 1 : length(profX.plat);
for iiiiA = 1:length(iiii)
  ii = iiii(iiiiA);
  jj = ii;

  %%% this has all beeen superseded by stuff below *****************
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

  %%% above has all been superseded by *******************************
  %%% this new code!!!!!!!!!!!!!!!!
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

  watercld = interp1(log10(profX.plevs(:,ii)),profX.clwc(:,ii),log10(plevs));
  icecld   = interp1(log10(profX.plevs(:,ii)),profX.ciwc(:,ii),log10(plevs));
  ptemp    = interp1(log10(profX.plevs(:,ii)),profX.ptemp(:,ii),log10(plevs));

  if ~exist('aa','var')
    aa = [];
  end
  aa = cloud_mean_press(aa,xcumsum,icecld,watercld,plevs,ii);

  %%% above has all been superseded by *******************************
  %%% this new code!!!!!!!!!!!!!!!!

  [wOUT,wT,wB,wPeak,wN,wmaxN,wminN] = boxshape(watercld,rGaussianCutoff); newwater = wOUT;
  [iOUT,iT,iB,iPeak,iN,imaxN,iminN] = boxshape(icecld,rGaussianCutoff);   newice   = iOUT;

  if iPrint > 0
    poink = [ii wN iN wPeak iPeak];
    %fprintf(1,'%5i %3i %3i | %8.6e %8.6e %8.6e | %8.6e %8.6e %8.6e \n',poink);
    fprintf(1,'%5i %3i %3i | %8.6e \n',poink);
  elseif iPrint < 0 & mod(jj,1000) == 0
    tnow = toc;
    fprintf(1,' processed %5i in %8.6f minutes\n',ii,tnow/60);
  end

  if iN > 3
    [iN,iOUT,iT,iB,iPeak] = combine_clouds4t3(iN,iOUT,iT,iB,iPeak,plevs);
  end
  if wN > 3
    [wN,wOUT,wT,wB,wPeak] = combine_clouds4t3(wN,wOUT,wT,wB,wPeak,plevs);
  end

  if iN > 2
    [iN,iOUT,iT,iB,iPeak] = combine_clouds3t2(iN,iOUT,iT,iB,iPeak,plevs);
    end
  if wN > 2
    [wN,wOUT,wT,wB,wPeak] = combine_clouds3t2(wN,wOUT,wT,wB,wPeak,plevs);
  end

  if ((iN == 1 & wN == 2) | (iN == 2 & wN == 1) | (iN == 2 & wN == 2))
    if iN == 2
      [iN,iOUT,iT,iB,iPeak] = combine_clouds2t1(iN,iOUT,iT,iB,iPeak,plevs);
    end
    if wN == 2
      [wN,wOUT,wT,wB,wPeak] = combine_clouds2t1(wN,wOUT,wT,wB,wPeak,plevs);  
    end
  end

  [cT,cB,cOUT,cngwat,cTYPE,iFound] = combine_clouds(...
              iN,iOUT,iT,iB,iPeak,wN,wOUT,wT,wB,wPeak,plevs,profX.plevs(:,ii));

  prof = put_into_prof(prof,profX,ii,jj,plevs,ptemp,iLevsVers,...
                       cT,cB,cOUT,cngwat,cTYPE,iFound);

  if iPrint > 0
    print_ecmwfcld2sartacld
  end

  if iDoPlot > 0
    plot_ecmwfcld2sartacld
  end

end    %% loop over iiiiA

%% put in the cloud cumulative fraction info, so that it can be used if necessary by "reset_cprtop"
prof.watercldX = aa.watercldX;
prof.icecldX   = aa.icecldX;
prof.watercldY = aa.watercldY;
prof.icecldY   = aa.icecldY;
