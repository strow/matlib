function plotclouds2(p1,p2,fig1,fig2)

%% almost like plotclouds but can handle two different input structire

%% ICE
figure(fig1); clf
oo1=find(p1.cfrac > 0 & p1.cngwat > 0 & p1.cprtop > 0 & p1.ctype == 201);
oo1=find(p1.cprtop > 0 & p1.ctype == 201);  
  plot(p1.rlat(oo1),p1.cprtop(oo1),'bo'); set(gca,'ydir','reverse')
hold on
oo2=find(p1.cfrac2 > 0 & p1.cngwat2 > 0 & p1.cprtop2 > 0 & p1.ctype2 == 201);
oo2=find(p1.cprtop2 > 0 & p1.ctype2 == 201); 
 plot(p1.rlat(oo2),p1.cprtop2(oo2),'bx'); set(gca,'ydir','reverse')
 xlabel('lat'); ylabel('ice cloud top (mb)');
fprintf(1,'I1 : %6i %6i \n',length(oo1),length(oo2))
oo1=find(p2.cfrac > 0 & p2.cngwat > 0 & p2.cprtop > 0 & p2.ctype == 201);
oo1=find(p2.cprtop > 0 & p2.ctype == 201);  
  plot(p2.rlat(oo1),p2.cprtop(oo1),'ro'); set(gca,'ydir','reverse')
hold on
oo2=find(p2.cfrac2 > 0 & p2.cngwat2 > 0 & p2.cprtop2 > 0 & p2.ctype2 == 201);
oo2=find(p2.cprtop2 > 0 & p2.ctype2 == 201); 
 plot(p2.rlat(oo2),p2.cprtop2(oo2),'rx'); set(gca,'ydir','reverse')
 xlabel('lat'); ylabel('ice cloud top (mb)');
fprintf(1,'I2 : %6i %6i \n',length(oo1),length(oo2))
hl = legend('p1ctype','p1ctype2','p2ctype','p2ctype2'); set(hl,'fontsize',10)

%% WATER
figure(fig2); clf
oo1=find(p1.cfrac > 0 & p1.cngwat > 0 & p1.cprtop > 0 & p1.ctype == 101);
oo1=find(p1.cprtop > 0 & p1.ctype == 101); 
  plot(p1.rlat(oo1),p1.cprtop(oo1),'bo'); set(gca,'ydir','reverse')
hold on
oo2=find(p1.cfrac2 > 0 & p1.cngwat2 > 0 & p1.cprtop2 > 0 & p1.ctype2 == 101);
oo2=find(p1.cprtop2 > 0 & p1.ctype2 == 101); 
  plot(p1.rlat(oo2),p1.cprtop2(oo2),'bx'); set(gca,'ydir','reverse')
  xlabel('lat'); ylabel('water cloud top (mb)');
fprintf(1,'W1 : %6i %6i \n',length(oo1),length(oo2))
oo1=find(p2.cfrac > 0 & p2.cngwat > 0 & p2.cprtop > 0 & p2.ctype == 101);
oo1=find(p2.cprtop > 0 & p2.ctype == 101); 
  plot(p2.rlat(oo1),p2.cprtop(oo1),'ro'); set(gca,'ydir','reverse')
hold on
oo2=find(p2.cfrac2 > 0 & p2.cngwat2 > 0 & p2.cprtop2 > 0 & p2.ctype2 == 101);
oo2=find(p2.cprtop2 > 0 & p2.ctype2 == 101); 
  plot(p2.rlat(oo2),p2.cprtop2(oo2),'rx'); set(gca,'ydir','reverse')
  xlabel('lat'); ylabel('water cloud top (mb)');
hold off
fprintf(1,'W2 : %6i %6i \n',length(oo1),length(oo2))
hl = legend('p1ctype','p1ctype2','p2ctype','p2ctype2'); set(hl,'fontsize',10)
