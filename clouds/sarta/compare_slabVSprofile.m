function y = compare_slabVSprofile(p,ii)
%% this script assume you have profile "p" and allows you to compare slab vs cloud profile for water and ice

% ii = input('Enter profile number : ');

nlev   = p.nlevs(ii);
p_prof = p.plevs(1:nlev,ii);
ice_prof   = p.ciwc(1:nlev,ii);  % in g/g
water_prof = p.clwc(1:nlev,ii);  % in g/g
cc_prof    = p.cc(1:nlev,ii);    % 0 < cc < 1

type1 = p.ctype(ii);
type2 = p.ctype2(ii);

if type1 == 101
  color1 = 'c';
elseif type1 == 201
  color1 = 'm';
end

if type2 == 101
  color2 = 'c';
elseif type2 == 201
  color2 = 'm';
end

norm = 1e6;

xice   = ice_prof;   xice   = cumsum(xice);
xwater = water_prof; xwater = cumsum(xwater);
plot(water_prof,p_prof,'b',...
     ice_prof,p_prof,'r',...
     xwater,p_prof,'b--',...
     xice,p_prof,'r--','linewidth',2);

plot(water_prof,p_prof,'b',...
     ice_prof,p_prof,'r','linewidth',2);

y.water_prof = water_prof;
y.ice_prof   = ice_prof;
y.plevs      = p_prof;

if p.cngwat(ii) > 0
  line([0 p.cngwat(ii)/norm],[p.cprtop(ii) p.cprtop(ii)],'color',color1)
  line([0 p.cngwat(ii)/norm],[p.cprbot(ii) p.cprbot(ii)],'color',color1)
  line([p.cngwat(ii)/norm p.cngwat(ii)/norm],[p.cprbot(ii) p.cprtop(ii)],'color',color1)
  shade2(gcf,0,p.cprbot(ii),p.cngwat(ii)/norm,p.cprtop(ii)-p.cprbot(ii),color1,0.25);
  if isfield(p,'icecldX')
    line([0 p.cngwat(ii)/norm],[p.icecldX(ii) p.icecldX(ii)])
  end
end

if p.cngwat2(ii) > 0
  line([0 p.cngwat2(ii)/norm],[p.cprtop2(ii) p.cprtop2(ii)],'color',color2)
  line([0 p.cngwat2(ii)/norm],[p.cprbot2(ii) p.cprbot2(ii)],'color',color2)
  line([p.cngwat2(ii)/norm p.cngwat2(ii)/norm],[p.cprbot2(ii) p.cprtop2(ii)],'color',color2)
  shade2(gcf,0,p.cprbot2(ii),p.cngwat2(ii)/norm,p.cprtop2(ii)-p.cprbot2(ii),color2,0.25);
  if isfield(p,'watercldX')
    line([0 p.cngwat2(ii)/norm],[p.watercldX(ii) p.watercldX(ii)])
  end
end

set(gca,'ydir','reverse');

%hl = legend('water','ice','waterslab','iceslab');
hl = legend('water','ice');
set(hl,'fontsize',12)

xlabel('mixing ratio g/g','fontsize',10); 
ylabel('pressure (mb)','fontsize',10)
%set(gca,'fontsize',10)

hx = axis;
axis([hx(1) hx(2) 100 1000]);

