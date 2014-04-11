function [site str] = transcom_match_plus(lat, lon)
% function [site str] = transcom_match_plus(lat, lon)
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
%   str{:}  - Cell array of strings with the list of regions (as bellow)
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
%
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
%   -------- Added Regions ------- 
%   23 - Arctic Lands (Greenland)
%   24 - Antartica
%   25 - Mediterranean/Arabian/Black/Caspian seas
%   26 - Great Lakes
%   
% See TranscomRegionMatrix.mat for the data array.
%   
% (see http://transcom.project.asu.edu/transcom03_protocol_basisMap.php)
% 
% Breno Imbiriba - 2013.10.17

% Call Base trasncom_match.m routine
site = transcom_match(lat, lon);

%% Add Arctic Lands as region 23
%site(site==0 & lat> 57) = 23;
%% Add Antartica as region 24
%site(site==0 & lat<-60) = 24;
%% Add Mediterranean region as region 25
%site(site==0 & lon>-20 & lon<70 & lat>0 & lat<60) = 25;
%% Add Great Lakes as region 26
%site(site==0 & lon>-150 & lon<-50) = 26;

% Now we have TransCom regions from 0 to 26 (Ignore 0 for now)

%    ----------- Lands -------------
  str{01} = 'North American Boreal';
  str{02} = 'North American Temperate';
  str{03} = 'South American Tropical';
  str{04} = 'South American Temperate';
  str{05} = 'Northern Africa';
  str{06} = 'Southern Africa';
  str{07} = 'Eurasian Boreal';
  str{08} = 'Eurasian Temperate';
  str{09} = 'Tropical Asia';
  str{10} = 'Australia';
  str{11} = 'Europe';
%         - '------- Oceans -----------
  str{12} = 'North Pacific Temperate';
  str{13} = 'West Pacific Tropics';
  str{14} = 'East Pacific Tropics';
  str{15} = 'South Pacific Temperate';
  str{16} = 'Northen Ocean';
  str{17} = 'North Atlantic Temperate';
  str{18} = 'Atlantic Tropics ';
  str{19} = 'South Atlantic Temperate';
  str{20} = 'Southern Ocean';
  str{21} = 'Indian Ocean';
  str{22} = 'South Indian Temperate';

%         - '--- Added Regions ------- 
  str{23} = 'Arctic Lands';
  str{24} = 'Antartica';
  str{25} = 'Mediterranean/Arabian/Black/Caspian seas';
  str{26} = 'Great Lakes';

end
