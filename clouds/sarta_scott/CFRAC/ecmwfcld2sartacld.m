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

if nlev > 20
  iSmooth = 2;
else
  iSmooth = 1;
end

global iDoPlot
iDoPlot = +1;   %% plot stuff .. ugh slow
iDoPlot = -1;   %% do not plot stuff

iPrint = -1;    %% do not print chirpy talky comments
iPrint = +1;    %% do not print chirpy talky comments

tic;

iStep = 1;
jj = 0;
for ii = 1 : iStep : length(profX.plat)
  jj = jj + 1;

  cfracW = sum(profX.clwc(:,ii).*profX.cc(:,ii))/sum(profX.clwc(:,ii)+1e-15);
  cfracI = sum(profX.ciwc(:,ii).*profX.cc(:,ii))/sum(profX.ciwc(:,ii)+1e-15);
  cfracIW = [cfracI cfracW];

  if nlev > 20
    [shiftedx,shiftedy,plevs]    =smooth1aa(1:nlev,profX.plevs(:,ii)',iSmooth);
    [shiftedx,shiftedy,watercld] =smooth1aa(1:nlev,profX.clwc(:,ii)' ,iSmooth);
    [shiftedx,shiftedy,icecld]   =smooth1aa(1:nlev,profX.ciwc(:,ii)' ,iSmooth);
    [shiftedx,shiftedy,cldfrac]  =smooth1aa(1:nlev,profX.cc(:,ii)' ,  iSmooth);
  else
    plevs    = profX.plevs(:,ii);
    watercld = profX.clwc(:,ii);
    icecld   = profX.ciwc(:,ii);
    cldfrac  = profX.cc(:,ii);
  end

  [wCC,wOUT,wT,wB,wPeak,wN,wmaxN,wminN] = boxshape(cldfrac,watercld); 
     newwater = wOUT;
  [iCC,iOUT,iT,iB,iPeak,iN,imaxN,iminN] = boxshape(cldfrac,icecld);   
     newice   = iOUT;

  if iPrint > 0
    poink = [ii wN iN wPeak iPeak];
    fprintf(1,'%5i %3i %3i | %8.6e %8.6e %8.6e | %8.6e %8.6e %8.6e \n',poink);
  elseif iPrint < 0 & mod(jj,1000) == 0
    tnow = toc;
    fprintf(1,' processed %5i in %8.6f minutes\n',ii,tnow/60);
  end

  if iN > 3
    [iCC,iN,iOUT,iT,iB,iPeak]=combine_clouds4t3(iCC,iN,iOUT,iT,iB,iPeak,plevs);
  end
  if wN > 3
    [wCC,wN,wOUT,wT,wB,wPeak]=combine_clouds4t3(wCC,wN,wOUT,wT,wB,wPeak,plevs);
  end

  if iN > 2
    [iCC,iN,iOUT,iT,iB,iPeak]=combine_clouds3t2(iCC,iN,iOUT,iT,iB,iPeak,plevs);
  end
  if wN > 2
    [wCC,wN,wOUT,wT,wB,wPeak]=combine_clouds3t2(wCC,wN,wOUT,wT,wB,wPeak,plevs);
  end

  if ((iN == 1 & wN == 2) | (iN == 2 & wN == 1) | (iN == 2 & wN == 2))
    if iN == 2
      [iCC,iN,iOUT,iT,iB,iPeak] = ...
         combine_clouds2t1(iCC,iN,iOUT,iT,iB,iPeak,plevs);
    end
    if wN == 2
      [wCC,wN,wOUT,wT,wB,wPeak] = ...
         combine_clouds2t1(wCC,wN,wOUT,wT,wB,wPeak,plevs);  
    end
  end

  [cFrac,cT,cB,cOUT,cngwat,cTYPE,iFound] = combine_clouds(...
                               iCC,iN,iOUT,iT,iB,iPeak,...
                               wCC,wN,wOUT,wT,wB,wPeak,...
                               plevs,profX.plevs(:,ii));

  %%%%%% put_into_pnew; nah!!!!!!
  put_into_prof;

  if iPrint > 0
    fprintf(1,' \n');
    if (length(cTYPE) >= 1)
      disp('type    toplev   botlev    cngwat*1000 kg/kg     cngwat g/m2     cfracSimple cfracNew');
      disp('----------------------------------------------------------------------------------------');
      for kk = 1 : length(cTYPE)
        if cTYPE(kk) == 'I'
          lala = [cT(kk),cB(kk),cngwat(kk)*1000,cc(kk),cfracIW(kk),cFrac(kk))];
          fprintf(1,'  I      %3i       %3i         %8.6f         %8.6f      %8.6f  %8.6f\n',lala)
      else
          lala = [cT(kk),cB(kk),cngwat(kk)*1000,cc(kk),cfracIW(kk),cFrac(kk))];
          fprintf(1,'  W      %3i       %3i         %8.6f         %8.6f      %8.6f  %8.6f\n',lala)
        end
      end
      disp('----------------------------------------------------------------------------------------');
    else
      fprintf(1,' no cld found \n');
    end
    fprintf(1,' \n');
  end

  if iDoPlot > 0
    figure(1);  clf
    plot(profX.ciwc(:,ii),profX.plevs(:,ii),'b',...
        icecld,plevs,'b--',newice,plevs,'bo-',...
        profX.clwc(:,ii),profX.plevs(:,ii),'r',watercld,plevs,'r--',...
                                      newwater,plevs,'ro-'); hold on
    plot(profX.ciwc(:,ii),profX.plevs(:,ii),'b',...
         profX.clwc(:,ii),profX.plevs(:,ii),'r',...
         cOUT,profX.plevs(:,ii),'k','LineWidth',3 )
    set(gca,'ydir','Reverse'); 
    grid
    hold off
    legend('Ice','SmoothIce','NewIce','Water','SmoothWater','NewWater',...
         'Location','NorthEast')
    pause(0.1);
  end

  subplot(121); plot(iCC,plevs,cldfrac,...
                     plevs,icecld/(max(icecld)+1e-10),plevs); 
    title('I');
  subplot(122); plot(wCC,plevs,cldfrac,...
                     plevs,watercld/(max(watercld)+1e-10),plevs); 
    title('W');
  pause

end

