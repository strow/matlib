jett = jet; jett(1,:) = 0.6;
ctmp = coast;
for dd = 20 : 24
  for gg = 1 : 240

    iYear = num2str(2007);
    iMonth = 02;
    if iMonth < 10
      iMonth = ['0' num2str(iMonth)];
    else
      iMonth = [    num2str(iMonth)];
      end
    iDay = dd;
    if iDay < 10
      iDay = ['0' num2str(iDay)];
    else
      iDay = [    num2str(iDay)];
      end
    iGran = gg;
    if iGran < 10
      iGran = ['00' num2str(iGran)];
    elseif iGran < 100
      iGran = ['0' num2str(iGran)];
    else
      iGran = [    num2str(iGran)];
      end

    fname = ['winds_' iYear '_' iMonth '_' iDay '_' iGran '.mat'];
    gname = ['winds ' iYear ' ' iMonth ' ' iDay ' ' iGran ];
    pname = ['winds_' iYear '_' iMonth '_' iDay '_' iGran];

    cname = ['/carrot/s1/sergio/SPECIALRTPFILES/' iYear '/' iMonth '/' iDay];
    cname = [cname '/resetcloudparam_optimum_A_Q_' iGran '_STreset_N3_DustOnly.mat'];

    ee = exist(fname,'file');
    cc = exist(cname,'file');
    if ee > 0 & cc > 0

      loader = ['load ' fname];
      eval(loader)

      figure(1)
      ZZ = U10m.^2 + V10m.^2; ZZ = sqrt(ZZ);
      ind = 1:1:120;ind = 1:5:120;
      scatter_coast(XX(:),YY(:),30,ZZ(:)); hold on
      quiver(XX(ind,ind),YY(ind,ind),U10m(ind,ind),V10m(ind,ind),'k'); 
      title([gname ' 10m spd']);
      hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
      axis([-15 45 -10 50]);
      caxis([0 20]); colorbar
      colormap(jett)

      loader = ['load ' cname];
      eval(loader)
      p = what2save;
      haha = find(p.raTau900 > 0.2);
      figure(2); scatter_coast(p.rlon(haha),p.rlat(haha),20,p.raTau900(haha));
      title('OD');

      XXX = XX(:); YYY = YY(:); ZZZ = ZZ(:);
      savespeed = [];
      for mm = 1 : length(haha)
        dist = (XXX-p.rlon(haha(mm))).^2 + (YYY-p.rlat(haha(mm))).^2; 
        dist = sqrt(dist);
        tada = find(dist <= 1);
        if length(tada) > 2
          xxsave(mm) = mean(XXX(tada));   yysave(mm) = mean(YYY(tada)); 
          savespeed     = [savespeed; ZZZ(tada)];    
          speedmean(mm) = mean(ZZZ(tada));
          speedstd(mm)  = std(ZZZ(tada));
          speedpeak(mm) = max(ZZZ(tada));
          end
        end
      figure(4); scatter_coast(xxsave,yysave,20,speedmean); title('Mean speed');

      array = [str2num(iGran)         mean(xxsave)           mean(yysave)          ...
               max(p.raTau900(haha))  mean(p.raTau900(haha)) std(p.raTau900(haha)) ...
               max(savespeed)         mean(savespeed)        std(savespeed)];

      fprintf(1,' %3s & %3i & %5.2f & %5.2f & %4.2f & %4.2f $\\pm$ %4.2f & %4.2f & %4.2f $\\pm$ %4.2f \n',iDay,array);
      %disp('ret to continue'); pause
      pause(1)
      end
    end
  end