%% for a post-processing version, see compare_slabVSprofile.m or compare_slabVSprofile2.m

figure(1); clf
lalap = profX.plevs(:,ii);
lalam = profX.ciwc(:,ii);
meanpI = sum(profX.plevs(:,ii).*profX.ciwc(:,ii))/(sum(profX.ciwc(:,ii))+1e-16);
plot(iOUT,plevs,'bx-',icecld,plevs,'c',profX.ciwc(:,ii),profX.plevs(:,ii),'k')
line([0 max(iOUT)*1.25],[prof.cprtop(ii) prof.cprtop(ii)],'color','k')
line([0 max(iOUT)*1.25],[prof.cprbot(ii) prof.cprbot(ii)],'color','k')
line([max(iOUT)*1.25 max(iOUT)*1.25],[prof.cprtop(ii) prof.cprbot(ii)],'color','k')
line([0 max(lalam)],[meanpI meanpI],'color','r','linewidth',2)
title('ICE')
set(gca,'ydir','reverse');

figure(2); clf
lalap = profX.plevs(:,ii);
lalam = profX.clwc(:,ii);
meanpI = sum(profX.plevs(:,ii).*profX.clwc(:,ii))/(sum(profX.clwc(:,ii))+1e-16);
plot(wOUT,plevs,'bx-',watercld,plevs,'c',profX.clwc(:,ii),profX.plevs(:,ii),'k')
line([0 max(wOUT)*1.25],[prof.udef(13,ii) prof.udef(13,ii)],'color','k')
line([0 max(wOUT)*1.25],[prof.udef(14,ii) prof.udef(14,ii)],'color','k')
line([max(wOUT)*1.25 max(wOUT)*1.25],[prof.udef(13,ii) prof.udef(14,ii)],'color','k')
line([0 max(lalam)],[meanpI meanpI],'color','r','linewidth',2)
title('WATER')
set(gca,'ydir','reverse');

whos cOUT

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
hl = legend('Ice','SmoothIce','NewIce','Water','SmoothWater','NewWater',...
     'Location','NorthEast');
set(hl,'fontsize',8);
title('BOTH')     

figure(4); clf
plot(profX.ptemp(:,ii),profX.plevs(:,ii)); title('T(z)');
set(gca,'ydir','Reverse'); 
xlim([180 320])
grid

for ijunki = 1 : 4
  figure(ijunki); ylim([min(profX.plevs(:,ii)) profX.spres(ii) + 10]);
end

figure(1); set(gca,'xscale','log')
figure(2); set(gca,'xscale','log')
figure(3); set(gca,'xscale','log')
figure(4); set(gca,'xscale','log')
figure(1); set(gca,'yscale','log')
figure(2); set(gca,'yscale','log')
figure(3); set(gca,'yscale','log')
figure(4); set(gca,'yscale','log')

pause(0.1);

