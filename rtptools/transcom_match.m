function site = transcom_match(lat, lon)
%function site = transcom_match(lat, lon)
%
% Using the TRANSCOM Terrestrial and Oceanic regions definitions, 
% at 1x1 degree resolution, assign a regions to each provided (lat,lon) pair.
% 
%
% INPUTS
%   lat(:),lon(:) - vectos of locations.
%
% OUTPUT
%   site(:) - TRANSCOM Region index of each location.
%
% TRANSCOM Regions:
%
%   ----------- Lands -------------
%    1 - North American Boreal
%    2 - North American Temperate
%    3 - South American Tropical
%    4 - South American Temperate
%    5 - Northern Africa
%    6 - Southern Africa
%    7 - Eurasian Boreal
%    8 - Eurasian Temperate
%    9 - Tropical Asia
%   10 - Australia
%   11 - Europe
%   ----------- Oceans -----------
%   12 - North Pacific Temperate
%   13 - West Pacific Tropics
%   14 - East Pacific Tropics
%   15 - South Pacific Temperate
%   16 - Northen Ocean
%   17 - North Atlantic Temperate
%   18 - Atlantic Tropics 
%   19 - South Atlantic Temperate
%   20 - Southern Ocean
%   21 - Indian Ocean
%   22 - South Indian Temperate
% 
%   -- Regions NOT Represented --
%    0 - Greenland
%    0 - Antartica
%    0 - Mediterranean/Arabian/Black/Caspian seas
%
% See TranscomRegionMatrix.mat for the data array.
% 
% (see http://transcom.project.asu.edu/transcom03_protocol_basisMap.php)
% 
% Paul Schou?

load TranscomRegionMatrix
ii = sub2ind(size(RegionMatrix),max(1,min(180,round(-lat+89.5))),mod(round(lon-0.5),360)+1);
site = RegionMatrix(ii);

