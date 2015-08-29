function [prof,profX] = ecmwfcld2sartacld(profIN,nlev,xcumsum,airslevels,airsheights);

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

iMakeIceWaterCld = -1; %% old style for taking eg cc/ciwc/clwc and smoothing
iMakeIceWaterCld = +1; %% new style for taking eg cc/ciwc/clwc and smoothing

tic;

jj = 0;
iiii = 1 : length(profX.plat);
for iiiiA = 1:length(iiii)
  ii = iiii(iiiiA);
  jj = ii;

  if iMakeIceWaterCld == -1
    old_style_smooth_cc_ciwc_clwc_to_water_ice_profile
  else
    new_style_smooth_cc_ciwc_clwc_to_water_ice_profile
  end
  
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

  cloud_combine_main_code
  
  prof = put_into_prof(prof,profX,ii,jj,plevs,ptemp,iLevsVers,...
                       cT,cB,cOUT,cngwat,cTYPE,iFound,airslevels,airsheights);

  if iPrint > 0
    print_ecmwfcld2sartacld
  end

  if iDoPlot > 0
    plot_ecmwfcld2sartacld
  end

%{
  %% debug
  if iiiiA == 1833
    iiiiA
    cT
    cB
    [iT iB wT wB]
    junky = iiiiA;
    [profIN.cprtop(junky) profIN.cprbot(junky) profIN.cprtop2(junky) profIN.cprbot2(junky) profIN.udef(13,junky) profIN.udef(14,junky)]
    [prof.cprtop(junky) prof.cprbot(junky) prof.cprtop2(junky) prof.cprbot2(junky) prof.udef(13,junky) prof.udef(14,junky)]
    plot(profX.plevs(:,ii),profX.clwc(:,ii),'bo-',plevs,watercld,'r-',plevs([wT wB]),[1 1]*wPeak,'kx-')
    keyboard
  end
%}

end    %% loop over iiiiA

%% put in the cloud cumulative fraction info, so that it can be used if necessary by "reset_cprtop"
prof.watercldX = aa.watercldX;
prof.icecldX   = aa.icecldX;
prof.watercldY = aa.watercldY;
prof.icecldY   = aa.icecldY;
