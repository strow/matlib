function [lat lon name pop] = worldcities(pop_range,lat_range,lon_range)
%function [lat lon name pop] = worldcities(pop_range,lat_range,lon_range)
%
% MaxMind's database of world cities and their populations
%

% To get the raw data download this file: 
%   http://www.maxmind.com/download/worldcities/worldcitiespop.txt.gz
% convert it to unicode using:
%   iconv -f ISO-8859-1 -t UTF-8 worldcitiespop.txt > worldcitiespop2.txt
% and read it in using:
% [Country,City,AccentCity,Region,Population,Latitude,Longitude]=textread('worldcitiespop2.txt','%s%s%s%s%d%f%f','delimiter',',','headerlines',1,'emptyvalue',NaN);
% s = ~isnan(Population);
% Latitude=single(Latitude(s));
% Longitude=single(Longitude(s));
% Population=single(Population(s));
% AccentCity=AccentCity(s);
% clear Country Region City


load worldcities.mat
s = ones(size(Latitude),'uint8');

if nargin > 2 & ~isempty(pop_range)
  if length(pop_range) == 1
    s = s & Population > pop_range;
  else
    s = s & Population > pop_range(1) & Population < pop_range(2);
  end
end
if nargin > 0 & ~isempty(lat_range)
  s = s & Latitude > lat_range(1) & Latitude < lat_range(2);
end
if nargin > 1 & ~isempty(lon_range)
  s = s & Longitude > lon_range(1) & Longitude < lon_range(2);
end
lat = Latitude(s);
lon = Longitude(s);
name = AccentCity(s);
pop = Population(s);

[name i]=unique(name);
lat = lat(i);
lon = lon(i);
pop = pop(i);
