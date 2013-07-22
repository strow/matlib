filename = '/carrot/s1/sergio/ECMWFJUNK/2007/02/';
filename = '/carrot/s1/sergio/ECMWFJUNK/2007/02/UAD02201200022012001';
filename = '/carrot/s1/sergio/ECMWFJUNK/2007/02/UAD02221200022212001';
filename = '/carrot/s1/sergio/ECMWFJUNK/2007/02/UAD02220000022203001';

minlat = 0;   maxlat = +60;
minlon = -20 + 180; maxlon = +40 + 180; 

minlat = -90;     maxlat = +90;
minlon = 0;       maxlon = +359;
gridsize = 0.50;

[head, profX] = ...
  readecmwf91_grid_winds(filename, minlat, maxlat, minlon, maxlon, gridsize);
profX.plon = npi2pi(profX.plon);

windplot;
