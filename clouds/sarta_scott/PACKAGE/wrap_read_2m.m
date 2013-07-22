xx = input('Enter [YY MM DD GG] : ');

strYear = num2str(xx(1));
iMonth = xx(2);
if iMonth < 10
  strMonth = ['0' num2str(iMonth)];
else
  strMonth = [    num2str(iMonth)];
  end
iDay = xx(3);
if iDay < 10
  strDay = ['0' num2str(iDay)];
else
  strDay = [    num2str(iDay)];
  end
iGran = xx(4);
if iGran < 10
  strGran = ['00' num2str(iGran)];
elseif iGran < 100
  strGran = ['0' num2str(iGran)];
else
  strGran = [    num2str(iGran)];
  end

rtpfile = ['/asl/data/rtprod/' strYear '/' strMonth '/'];
rtpfile = [rtpfile strDay '/allfov' strGran '.rtp'];
ee = exist(rtpfile,'file');
if ee == 0
  error('rtpfile does not exist')
else
  [h,ha,p,pa] = rtpread(rtpfile);
  end

xxstr = ...
  [num2str(xx(1)) ' ' num2str(xx(2)) ' ' num2str(xx(3)) ' ' num2str(xx(4))];
homedir = pwd;
cd /home/sergio/Bin
runbin = ['!find_granule_ecmwf91 ' xxstr ' >& ugh'];
eval(runbin)
sedder = ['!sed -e ''1,5d'' -e ''7,7d'' ugh > ugh1'];
eval(sedder); 
fid = fopen('ugh1');
tline = fgetl(fid);
disp(tline)
fclose(fid);
rmer = ['!/bin/rm ugh1 ugh']; eval(rmer);
cder = ['cd ' homedir]; eval(cder);

filename = tline;

minlat = 0;   maxlat = +60;
minlon = -20 + 180; maxlon = +40 + 180; 

minlat = -90;     maxlat = +90;
minlon = 0;       maxlon = +359.5;
gridsize = 0.5;

%%  find(profX.plon >= -15 & profX.plon < 45 & ...
%%       profX.plat >= -10 & profX.plat < 50);
minlat = -20;     maxlat = +50;
minlon = -20;     maxlon = +50;  
gridsize = 0.25;
minlon = zero22pi(minlon);
maxlon = zero22pi(maxlon);

minlat = -20;     maxlat = +50;
minlon = 0;       maxlon = +359.75;
gridsize = 0.25;

[head, profX] = ...
  readecmwf91_grid_winds(filename, minlat, maxlat, minlon, maxlon, gridsize);
profX.plon = npi2pi(profX.plon);

ii = find(profX.plon >= min(p.rlon)-2 & profX.plon < max(p.rlon)+2 & ...
          profX.plat >= min(p.rlat)-2 & profX.plat < max(p.rlat)+2);
fprintf(1,'%7i to be saved out of %7i \n',length(ii),length(profX.stemp));
plot(profX.plon(ii),profX.plat(ii),'.')
profX.plat = profX.plat(ii);
profX.plon = profX.plon(ii);
profX.T2 = profX.T2(ii);
profX.U10 = profX.U10(ii);
profX.V10 = profX.V10(ii);
profX.stemp = profX.stemp(ii);
profX.spres = profX.spres(ii);
profX.nlevs = profX.nlevs(ii);
profX.cfrac = profX.cfrac(ii);
profX.plevs = profX.plevs(:,ii);
profX.ptemp = profX.ptemp(:,ii);
profX.cc = profX.cc(:,ii);
profX.uspd = profX.uspd(:,ii);
profX.vspd = profX.vspd(:,ii);

for ii = 1 : length(p.stemp)
  if mod(ii,1000) == 0
    fprintf(1,'%5i out of %5i \n',ii,length(p.stemp));
    end
  dist = (p.rlon(ii) - profX.plon).^2 + (p.rlat(ii) - profX.plat).^2;
  dist = sqrt(dist);
  jj = find(dist == min(dist)); jj = jj(1);
  p2m.plevs(:,ii) = profX.plevs(:,jj);
  p2m.ptemp(:,ii) = profX.ptemp(:,jj);
  p2m.uspd(:,ii) = profX.uspd(:,jj);
  p2m.vspd(:,ii) = profX.vspd(:,jj);
  p2m.T2(ii)  = profX.T2(jj);
  p2m.U10(ii) = profX.U10(jj);
  p2m.V10(ii) = profX.V10(jj);
  p2m.stemp(ii) = profX.stemp(jj);
  p2m.spres(ii) = profX.spres(jj);
  p2m.rlon(ii) = p.rlon(ii);
  p2m.rlat(ii) = p.rlat(ii);
  end

saver = ['save t2m_winds_' strYear '_' strMonth '_' strDay '_' strGran ' p2m'];
eval(saver)

%windplot;
%windplot_wun;
