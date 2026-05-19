if iN > 3
  [iN,iOUT,iT,iB,iPeak] = combine_clouds4t3(iN,iOUT,iT,iB,iPeak,plevs,airslevels,airslayers,airsheights);
end
if wN > 3
  [wN,wOUT,wT,wB,wPeak] = combine_clouds4t3(wN,wOUT,wT,wB,wPeak,plevs,airslevels,airslayers,airsheights);
end

if iN > 2
  [iN,iOUT,iT,iB,iPeak] = combine_clouds3t2(iN,iOUT,iT,iB,iPeak,plevs,airslevels,airslayers,airsheights);
end
if wN > 2
  [wN,wOUT,wT,wB,wPeak] = combine_clouds3t2(wN,wOUT,wT,wB,wPeak,plevs,airslevels,airslayers,airsheights);
end

if ((iN == 1 & wN == 2) | (iN == 2 & wN == 1) | (iN == 2 & wN == 2))
  if iN == 2
    [iN,iOUT,iT,iB,iPeak] = combine_clouds2t1(iN,iOUT,iT,iB,iPeak,plevs,airslevels,airslayers,airsheights);
  end
  if wN == 2
    [wN,wOUT,wT,wB,wPeak] = combine_clouds2t1(wN,wOUT,wT,wB,wPeak,plevs,airslevels,airslayers,airsheights);  
  end
end

[cT,cB,cOUT,cngwat,cTYPE,iFound] = combine_clouds(...
            iN,iOUT,iT,iB,iPeak,wN,wOUT,wT,wB,wPeak,plevs,profX.plevs(:,ii),airslevels,airslayers,airsheights);

iTT = iT; if length(iTT) == 0; iTT = -1; end
iBB = iB; if length(iBB) == 0; iBB = -1; end
wTT = wT; if length(wTT) == 0; wTT = -1; end
wBB = wB; if length(wBB) == 0; wBB = -1; end
%fprintf(1,'>>>> here 1 ii iT wT cT iB wB cB = %5i %4i %4i %4i %4i %4i %4i \n',ii,iTT,wTT,cT(1),iBB,wBB,cB(1))

%% p440 and cut440 come from new_style_smooth_cc_ciwc_clwc_to_water_ice_profile
if p440 < 1013
  %% if ice_water_separator == 0 then p440 = 9999 since we do not want any changing
  if (length(cTYPE) >= 1)
    for kk = 1 : length(cTYPE)
      if cTYPE(kk) == 'I' & plevs(cT(kk)) > p440
        cT(kk) = cut440 - 10;
        cT(kk) = cut440 - 1;	
        if iPrint > 0
          fprintf(1,'warning had to reset ICE cloudtop, to make it less than %8.3f mb \n',p440)
        end
      elseif cTYPE(kk) == 'W' & plevs(cT(kk)) < p440
        cT(kk) = cut440 + 1;
        %cB(kk) = cB(kk)+3;
        if cB(kk) < cT(kk)
          cB(kk) = cT(kk) + 3;
        end
        if iPrint > 0
          fprintf(1,'warning had to reset WATER cloudtop, to make it more than %8.3f mb \n',p440)      
        end
      end
    end
  end
end

iTT = iT; if length(iTT) == 0; iTT = -1; end
iBB = iB; if length(iBB) == 0; iBB = -1; end
wTT = wT; if length(wTT) == 0; wTT = -1; end
wBB = wB; if length(wBB) == 0; wBB = -1; end
%fprintf(1,'>>>> here 2 ii iT wT cT iB wB cB = %5i %4i %4i %4i %4i %4i %4i \n',ii,iTT,wTT,cT(1),iBB,wBB,cB(1))
