function rcalcdiff = compare_profstructs(h,p1,p2,iIndex);

%% simple routine to compare profile "iIndex" of input rtp structs  p1 and p2

rcalcdiff = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
xyz = [p1.cfrac(iIndex) p1.cpsize(iIndex) p1.cngwat(iIndex) p1.cprtop(iIndex) p1.cprbot(iIndex)];
fprintf(1,'prof1, cld1 : cfrac cpsze cngwat cprtop cprbot = %8.6f %8.6f %8.6f %8.6f %8.6f \n',xyz);
hh=line([0 p1.cngwat(iIndex)],[p1.cprtop(iIndex) p1.cprtop(iIndex)],'color','b'); 
  set(hh,'LineStyle','-'); set(hh,'Marker','o');
hh=line([0 p1.cngwat(iIndex)],[p1.cprbot(iIndex) p1.cprbot(iIndex)],'color','b'); 
  set(hh,'LineStyle','-'); set(hh,'Marker','o');
hh=line([p1.cngwat(iIndex) p1.cngwat(iIndex)],[p1.cprbot(iIndex) p1.cprtop(iIndex)],'color','b'); 
  set(hh,'LineStyle','-'); set(hh,'Marker','o');

xyz = [p1.cfrac2(iIndex) p1.cpsize2(iIndex) p1.cngwat2(iIndex) p1.cprtop2(iIndex) p1.cprbot2(iIndex)];
fprintf(1,'prof1, cld2 : cfrac cpsze cngwat cprtop cprbot = %8.6f %8.6f %8.6f %8.6f %8.6f \n',xyz);
hh=line([0 p1.cngwat2(iIndex)],[p1.cprtop2(iIndex) p1.cprtop2(iIndex)],'color','r');
  set(hh,'LineStyle','-'); set(hh,'Marker','o');
hh=line([0 p1.cngwat2(iIndex)],[p1.cprbot2(iIndex) p1.cprbot2(iIndex)],'color','r');
  set(hh,'LineStyle','-'); set(hh,'Marker','o');
hh=line([p1.cngwat2(iIndex) p1.cngwat2(iIndex)],[p1.cprbot2(iIndex) p1.cprtop2(iIndex)],'color','r');
  set(hh,'LineStyle','-'); set(hh,'Marker','o');

set(gca,'ydir','reverse')
%%%%%%%%%%%%%%%%%%%%%%%%%
xyz = [p2.cfrac(iIndex) p2.cpsize(iIndex) p2.cngwat(iIndex) p2.cprtop(iIndex) p2.cprbot(iIndex)];
fprintf(1,'prof2, cld1 : cfrac cpsze cngwat cprtop cprbot = %8.6f %8.6f %8.6f %8.6f %8.6f \n',xyz);
line([0 p2.cngwat(iIndex)],[p2.cprtop(iIndex) p2.cprtop(iIndex)],'color','c');
line([0 p2.cngwat(iIndex)],[p2.cprbot(iIndex) p2.cprbot(iIndex)],'color','c');
line([p2.cngwat(iIndex) p2.cngwat(iIndex)],[p2.cprbot(iIndex) p2.cprtop(iIndex)],'color','c');

xyz = [p2.cfrac2(iIndex) p2.cpsize2(iIndex) p2.cngwat2(iIndex) p2.cprtop2(iIndex) p2.cprbot2(iIndex)];
fprintf(1,'prof2, cld2 : cfrac cpsze cngwat cprtop cprbot = %8.6f %8.6f %8.6f %8.6f %8.6f \n',xyz);
line([0 p2.cngwat2(iIndex)],[p2.cprtop2(iIndex) p2.cprtop2(iIndex)],'color','m');
line([0 p2.cngwat2(iIndex)],[p2.cprbot2(iIndex) p2.cprbot2(iIndex)],'color','m');
line([p2.cngwat2(iIndex) p2.cngwat2(iIndex)],[p2.cprbot2(iIndex) p2.cprtop2(iIndex)],'color','m');

boo = max([p1.cngwat(iIndex) p1.cngwat2(iIndex) p2.cngwat(iIndex) p2.cngwat2(iIndex)]);
if (boo > eps)
  axis([0 boo*1.1 0 1000]);
end

disp(' ')
boo = [p1.ctype(iIndex) p1.ctype2(iIndex) p1.cfrac12(iIndex)];
fprintf(1,'prof1 : ctype ctype2 cfrac12 = %8.6f %8.6f %8.6f \n',boo);
boo = [p2.ctype(iIndex) p2.ctype2(iIndex) p2.cfrac12(iIndex)];
fprintf(1,'prof2 : ctype ctype2 cfrac12 = %8.6f %8.6f %8.6f \n',boo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
boo = [p1.spres(iIndex) p1.stemp(iIndex) p1.scanang(iIndex) p1.solzen(iIndex) p1.rlat(iIndex) p1.rlon(iIndex)];
fprintf(1,'prof1 : spres stemp scanang solzen lat lon = %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',boo);
boo = [p2.spres(iIndex) p2.stemp(iIndex) p2.scanang(iIndex) p2.solzen(iIndex) p2.rlat(iIndex) p2.rlon(iIndex)];
fprintf(1,'prof2 : spres stemp scanang solzen lat lon = %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',boo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlevs1 = double(p1.nlevs(iIndex));
nlevs2 = double(p2.nlevs(iIndex));


figure(2); clf
subplot(121); 
  plot(p1.ptemp(1:nlevs1,iIndex),1:nlevs1,'bo-',p2.ptemp(1:nlevs2,iIndex),1:nlevs2,'rx-')
  set(gca,'ydir','reverse'); title('T(z)')
subplot(122); 
  semilogx(p1.gas_1(1:nlevs1,iIndex),1:nlevs1,'bo-',p2.gas_1(1:nlevs2,iIndex),1:nlevs2,'rx-')
  set(gca,'ydir','reverse'); title('WV(z)')
%[nansum(p1.ptemp(1:nlevs1,iIndex) - p2.ptemp(1:nlevs2,iIndex)) ...
% (nansum(p1.gas_1(1:nlevs1,iIndex)./p2.gas_1(1:nlevs2,iIndex)))/nlevs1]

figure(3); clf
subplot(131); 
  plot(p1.ciwc(1:nlevs1,iIndex),p1.plevs(1:nlevs1,iIndex),'bo-',...
       p2.ciwc(1:nlevs2,iIndex),p2.plevs(1:nlevs2,iIndex),'rx-')
  set(gca,'ydir','reverse'); title('CIWC(z)')
subplot(132); 
  plot(p1.clwc(1:nlevs1,iIndex),p1.plevs(1:nlevs1,iIndex),'bo-',...,
       p2.clwc(1:nlevs2,iIndex),p2.plevs(1:nlevs2,iIndex),'rx-')
  set(gca,'ydir','reverse'); title('CLWC(z)')
subplot(133); 
  plot(p1.cc(1:nlevs1,iIndex),p1.plevs(1:nlevs1,iIndex),'bo-',...,
       p2.cc(1:nlevs2,iIndex),p2.plevs(1:nlevs2,iIndex),'rx-')
  set(gca,'ydir','reverse'); title('CC(z)')
%[(nansum(p1.ciwc(1:nlevs1,iIndex)./p2.ciwc(1:nlevs2,iIndex)))/nlevs1 ...
% (nansum(p1.clwc(1:nlevs1,iIndex)./p2.clwc(1:nlevs2,iIndex)))/nlevs1 ...
%  nansum(p1.cc(1:nlevs1,iIndex) - p2.cc(1:nlevs2,iIndex))]

%figure(5)
%  plot(p1.clwc(1:nlevs1,iIndex)./p2.clwc(1:nlevs2,iIndex),p1.plevs(1:nlevs1,iIndex),'bo-')
%  set(gca,'ydir','reverse'); title('CLWC(z)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf
lala = ones(size(h.vchan)) * nan;

if isfield(p1,'robs1')
  robs_prof1 = p1.robs1(:,iIndex);
else
  robs_prof1 = lala;
end

if isfield(p2,'robs1')
  robs_prof2 = p2.robs1(:,iIndex);
else
  robs_prof2 = lala;
end

if isfield(p1,'rcalc')
  rcal_prof1 = p1.rcalc(:,iIndex);
else
  rcal_prof1 = lala;
end

if isfield(p2,'rcalc')
  rcal_prof2 = p2.rcalc(:,iIndex);
else
  rcal_prof2 = lala;
end

if isfield(p1,'sarta_rclearcalc')
  rclr_prof1 = p1.sarta_rclearcalc(:,iIndex);
else
  rclr_prof1 = lala;
end
if isfield(p2,'sarta_rclearcalc')
  rclr_prof2 = p2.sarta_rclearcalc(:,iIndex);
else
  rclr_prof2 = lala;
end

rall1 = [robs_prof1 rcal_prof1 rclr_prof1];
rall2 = [robs_prof2 rcal_prof2 rclr_prof2];
tall1 = rad2bt(h.vchan,rall1);
tall2 = rad2bt(h.vchan,rall2);
plot(h.vchan,tall1); hold on
plot(h.vchan,tall2,'--'); hold off
axis([min(h.vchan)-5 max(h.vchan)+5 180 300]); grid
title('BGR = obs cld clr; solid/dash = prof1/2')

woo = find(h.ichan == 1291);
if length(woo) == 1
  boo1 = tall1(woo,:);
  fprintf(1,'prof1 BT1231     obs  cld  clr = %8.6f %8.6f %8.6f \n',boo1);
  boo2 = tall2(woo,:);
  fprintf(1,'prof2 BT1231     obs  cld  clr = %8.6f %8.6f %8.6f \n',boo2);
  fprintf(1,'DBT prof1-prof2  obs  cld  clr = %8.6f %8.6f %8.6f \n',boo1-boo2);
else
  disp('no 1231 cm-1 channel');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
  disp(' ')
  boo = [nansum(p1.sarta_lvlODice(:,iIndex)) nansum(p1.sarta_lvlODwater(:,iIndex)) p1.sarta_lvl_iceOD_1(iIndex) p1.sarta_lvl_waterOD_1(iIndex)];
  fprintf(1,'prof1 aux : lvlODice lvlODwater lvlODice1 lvlODwater1 = %8.6f %8.6f %8.6f %8.6f \n',boo);

  boo = [nansum(p2.sarta_lvlODice(:,iIndex)) nansum(p2.sarta_lvlODwater(:,iIndex)) p2.sarta_lvl_iceOD_1(iIndex) p2.sarta_lvl_waterOD_1(iIndex)];
  fprintf(1,'prof2 aux : lvlODice lvlODwater lvlODice1 lvlODwater1 = %8.6f %8.6f %8.6f %8.6f \n',boo);
%}