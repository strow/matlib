if iN > 3
  [iN,iOUT,iT,iB,iPeak] = combine_clouds4t3(iN,iOUT,iT,iB,iPeak,plevs,airslevels,airsheights);
end
if wN > 3
  [wN,wOUT,wT,wB,wPeak] = combine_clouds4t3(wN,wOUT,wT,wB,wPeak,plevs,airslevels,airsheights);
end

if iN > 2
  [iN,iOUT,iT,iB,iPeak] = combine_clouds3t2(iN,iOUT,iT,iB,iPeak,plevs,airslevels,airsheights);
end
if wN > 2
  [wN,wOUT,wT,wB,wPeak] = combine_clouds3t2(wN,wOUT,wT,wB,wPeak,plevs,airslevels,airsheights);
end

if ((iN == 1 & wN == 2) | (iN == 2 & wN == 1) | (iN == 2 & wN == 2))
  if iN == 2
    [iN,iOUT,iT,iB,iPeak] = combine_clouds2t1(iN,iOUT,iT,iB,iPeak,plevs,airslevels,airsheights);
  end
  if wN == 2
    [wN,wOUT,wT,wB,wPeak] = combine_clouds2t1(wN,wOUT,wT,wB,wPeak,plevs,airslevels,airsheights);  
  end
end

[cT,cB,cOUT,cngwat,cTYPE,iFound] = combine_clouds(...
            iN,iOUT,iT,iB,iPeak,wN,wOUT,wT,wB,wPeak,plevs,profX.plevs(:,ii),airslevels,airsheights);

if (length(cTYPE) >= 1)
  for kk = 1 : length(cTYPE)
    if cTYPE(kk) == 'I' & plevs(cT(kk)) > 440
      cT(kk) = cut440 - 10;
      if iPrint > 0
        disp('warning had to reset ICE cloudtop, to make it less than 440 mb');
      end
    elseif cTYPE(kk) == 'W' & plevs(cT(kk)) < 440
      cT(kk) = cut440 + 1;
      cB(kk) = cB(kk)+3;
      if iPrint > 0	
        disp('warning had to reset WATER cloudtop, to make it more than 440 mb');
      end
    end
  end
end
