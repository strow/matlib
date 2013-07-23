xx = input('Enter [YY MM DD GG] : ');

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

%%  find(profX.plon >= -15 & profX.plon < 45 & profX.plat >= -10 & profX.plat < 50);
minlat = -20;     maxlat = +50;
minlon = -20;     maxlon = +50;  
gridsize = 0.25;
minlon = zero22pi(minlon);
maxlon = zero22pi(maxlon);

minlat = -20;     maxlat = +50;
minlon = 0;       maxlon = +359.75;
gridsize = 0.25;

[head, profX] = ...
readecmwf91_grid_winds_UVW(filename, minlat, maxlat, minlon, maxlon, gridsize);
profX.plon = npi2pi(profX.plon);

error('haha');

%windplot;
windplot_wun_UVW;
