delp = 0:10:1000;
delp = 0:50:1000;

figure(1); clf
oo11 = find(prof.ctype == 201 & prof.cngwat > 0);
%  plot(delp,hist(prof.cprtop(oo11),delp),'bo-'); grid;
oo12 = find(prof.ctype2 == 201 & prof.cngwat2 > 0);
%  plot(delp,hist(prof.cprtop2(oo12),delp),'bo-'); grid

oo21 = find(prof.ctype == 101 & prof.cngwat > 0);
%  plot(delp,hist(prof.cprtop(oo21),delp),'ro-'); grid  
oo22 = find(prof.ctype2 == 101 & prof.cngwat2 > 0);
%  plot(delp,hist(prof.cprtop2(oo22),delp),'o-'); grid

plot(delp,hist(prof.cprtop(oo11),delp),'bo-',delp,hist(prof.cprtop2(oo12),delp),'cx-',...
     delp,hist(prof.cprtop(oo21),delp),'ro-',delp,hist(prof.cprtop2(oo22),delp),'mx-','linewidth',2); grid
hl=legend('ice1','ice2','water1','water2'); set(hl,'fontsize',10);

%plot(prof.ctype,prof.cprtop,'.')
%plot(prof.ctype2,prof.cprtop2,'.')
%plot(prof.ciwc,prof.plevs)
