figure(1); clf
lalap = profX.plevs(:,ii);
lalam = profX.ciwc(:,ii);
meanpI = sum(profX.plevs(:,ii).*profX.ciwc(:,ii))/(sum(profX.ciwc(:,ii))+1e-16);
plot(plevs,iOUT,plevs,icecld,profX.plevs(:,ii),profX.ciwc(:,ii),'o-')
line([prof.cprtop(ii) prof.cprtop(ii)],[0 max(iOUT)*1.25],'color','k')
line([prof.cprbot(ii) prof.cprbot(ii)],[0 max(iOUT)*1.25],'color','k')
line([prof.cprtop(ii) prof.cprbot(ii)],[max(iOUT)*1.25 max(iOUT)*1.25],...
    'color','k')
line([meanpI meanpI],[0 max(lalam)],'color','r','linewidth',2)

figure(2); clf
lalap = profX.plevs(:,ii);
lalam = profX.clwc(:,ii);
meanpI = sum(profX.plevs(:,ii).*profX.clwc(:,ii))/(sum(profX.clwc(:,ii))+1e-16);
plot(plevs,wOUT,plevs,watercld,profX.plevs(:,ii),profX.clwc(:,ii),'o-')
line([prof.udef(13,ii) prof.udef(13,ii)],[0 max(wOUT)*1.25],'color','k')
line([prof.udef(14,ii) prof.udef(14,ii)],[0 max(wOUT)*1.25],'color','k')
line([prof.udef(13,ii) prof.udef(14,ii)],[max(wOUT)*1.25 max(wOUT)*1.25],...
    'color','k')
line([meanpI meanpI],[0 max(lalam)],'color','r','linewidth',2)


figure(3);  clf
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
