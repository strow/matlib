function plotclouds(p1,fig1,fig2,thestr)

%% ICE
figure(fig1); clf
oo1=find(p1.cfrac > 0 & p1.cngwat > 0 & p1.cprtop > 0 & p1.ctype == 201);  
  plot(p1.rlat(oo1),p1.cprtop(oo1),'bo'); set(gca,'ydir','reverse')
hold on
oo2=find(p1.cfrac2 > 0 & p1.cngwat2 > 0 & p1.cprtop2 > 0 & p1.ctype2 == 201); 
 plot(p1.rlat(oo2),p1.cprtop2(oo2),'cx'); set(gca,'ydir','reverse')
 xlabel('lat'); ylabel('ice cloud top (mb)');
hold off
fprintf(1,'I : %6i %6i \n',length(oo1),length(oo2))
if nargin == 4
  title(thestr);
end
axis([-90 +90 0 1000])

%% WATER
figure(fig2); clf
oo1=find(p1.cfrac > 0 & p1.cngwat > 0 & p1.cprtop > 0 & p1.ctype == 101); 
  plot(p1.rlat(oo1),p1.cprtop(oo1),'ro'); set(gca,'ydir','reverse')
hold on
oo2=find(p1.cfrac2 > 0 & p1.cngwat2 > 0 & p1.cprtop2 > 0 & p1.ctype2 == 101); 
  plot(p1.rlat(oo2),p1.cprtop2(oo2),'mx'); set(gca,'ydir','reverse')
  xlabel('lat'); ylabel('water cloud top (mb)');
hold off
fprintf(1,'W : %6i %6i \n',length(oo1),length(oo2))
if nargin == 4
  title(thestr);
end
axis([-90 +90 0 1000])