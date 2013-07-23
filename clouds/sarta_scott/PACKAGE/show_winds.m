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
    ee = exist(fname,'file');
    if ee > 0
      loader = ['load ' fname];
      eval(loader)

      figure(1)
      ZZ = U0.^2 + V0.^2; ZZ = sqrt(ZZ);
      ind = 1:1:120;ind = 1:5:120;
      scatter_coast(XX(:),YY(:),30,ZZ(:)); hold on
      quiver(XX(ind,ind),YY(ind,ind),U0(ind,ind),V0(ind,ind),'k'); 
      title([gname ' gnd spd']);
      hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
      axis([-15 45 -10 50]);
      caxis([0 20]); colorbar
      colormap(jett)

      figure(2)
      ZZ = U10m.^2 + V10m.^2; ZZ = sqrt(ZZ);
      ind = 1:1:120;ind = 1:5:120;
      scatter_coast(XX(:),YY(:),30,ZZ(:)); hold on
      quiver(XX(ind,ind),YY(ind,ind),U10m(ind,ind),V10m(ind,ind),'k'); 
      title([gname ' 10m spd']);
      hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
      axis([-15 45 -10 50]);
      caxis([0 20]); colorbar
      colormap(jett)

      figure(3)
      ZZ = U3km.^2 + V3km.^2; ZZ = sqrt(ZZ);
      ind = 1:1:120;ind = 1:5:120;
      scatter_coast(XX(:),YY(:),30,ZZ(:)); hold on
      quiver(XX(ind,ind),YY(ind,ind),U3km(ind,ind),V3km(ind,ind),'k'); 
      title([gname ' 2.5 km spd']);
      hold on; plot(ctmp(:,2),ctmp(:,1),'k','LineWIdth',2); hold off
      axis([-15 45 -10 50]);
      caxis([0 40]); colorbar
      colormap(jett)
     
      iPrint = +1;
      if iPrint > 0
        printer = ['figure(2); print -dpng ' pname]
        eval(printer);
        end
      %disp('ret to continue'); pause
      pause(2)
      end
    end
  end