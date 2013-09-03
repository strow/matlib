function [isiteind, isitenum] = fixedsite(lat, lon, range_km);

% function [isiteind, isitenum] = fixedsite(lat, lon, range_km);
%
% Find lat/lon within range of fixed sites.
%
% Input:
%    lat      : [1 x n] latitude -90 to 90
%    lon      : [1 x n] longitude -180 to 360
%    range_km : [1 x 1] max range (km)
%
% Output:
%    isiteind : [1 x m] indices of lat/lon within range
%    isitenum : [1 x m] site number
%

% Created: 11 April 2007, Scott Hannon
% Update: 20 April 2007, S.Hannon - fix "dist" conversion to km (was mm)
% Update:  3 March 2011, Paul Schou - added extra sites
%            2013.05.08, B.I. - removed extra sites (LLS request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Fixed sites used by AIRS
SiteLatLon=[ ...
    27.12   26.10; %  1 = Egypt 1
   -24.50  137.00; %  2 = Simpson Desert
   -75.10  123.40; %  3 = Dome Concordia, 3200 m elevation
     1.50  290.50; %  4 = Mitu, Columbia / Brazil tropical forest
     3.50   14.50; %  5 = Boumba, S.E. Cameroon
    38.50  244.30; %  6 = Railroad Valley, NV
    36.60  262.50; %  7 = ARM-SGP (southern great plains), OK
    -2.00  147.40; %  8 = Manus, Bismarck Archipelago
    -0.50  166.60; %  9 = ARM-TWP (tropical western pacific) Nauru, Micronesia
    90.00    0.00; % 10 = north pole
   -90.00    0.00; % 11 = south pole
    61.15   73.37; % 12 = Siberian tundra (Surgut)
    23.90  100.50; % 13 = Hunnan rain forest
    71.32  203.34; % 14 = Barrow, Alaska/ARM-NSA (north slope alaska)
    70.32  203.33; % 15 = Atqusuk, Alaska
   -12.42  130.89; % 16 = Darwin, Australia
    36.75  100.33; % 17 = Lake Qinhai
    40.17   94.33; % 18 = Dunhuang, Gobi desert
   -15.88  290.67; % 19 = Lake Titicaca
    39.10  239.96; % 20 = Lake Tahoe, CA
    31.05   57.65];% 21 = LUT Desert

% include global view sites
% load gvsites.mat
% SiteLatLon(101:100+length(gv_lat),:)=[gv_lat gv_lon];


% Convert Longitudes from 0:360 into -180:180
i360=find(SiteLatLon(:,2)>180);
SiteLatLon(i360,2) = SiteLatLon(i360,2)-360;

% If no input arguments, return lats and lons
if nargin == 0
  isiteind = SiteLatLon(:,1);
  isitenum = SiteLatLon(:,2);
  return
end

nSite = length(SiteLatLon);

% Approximately 111 km per 1 degree latitude
range_deg = range_km/111;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make up temporary work arrays
nin = length(lat);
dlatmax = range_deg; % latitude range

% Locate invalid latitudes 
igood = find(~isnan(lat));

% Clear up work arrays
sind = zeros(1,nin);
snum = zeros(1,nin);

% Loop over the fixed sites
for isite=1:nSite

   % skip over sites that are not assigned
   if isequal(SiteLatLon(isite,:),[0 0])
      continue;
   end
   
   slat = SiteLatLon(isite,1);
   slon = SiteLatLon(isite,2);

   % Use matlab "distance" command
   isel = distance(slat,slon,lat(igood),lon(igood)) < range_deg;

   % Mark positive hits
   sind(igood(isel)) = 1;
   snum(igood(isel)) = isite;
  
   % Remove them from the list - Show only the first match
   igood(isel) = [];  % Surprise syntax for me!

end % for nSite

% Assign output arrays
isiteind = find(sind == 1);
isitenum = snum(isiteind);

%%% end of function %%%
