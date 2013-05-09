function [efreq,emis] = emis_DanZhou(lat,lon,rtime,taiyear);

% function [efreq,emis] = emis_DanZhou(lat,lon,rtime,taiyear);
%
% Read Dan Zhous IASI emissivity database for 2007/07 to 2008/06
% and interpolate it for time of year and reduced the frequency
% point spacing (eg for use with RTP).  The year is ignored, only
% the time of year is used in the temporal interpolation.
%
% Input:
%    lat = [1 x nobs] latitude {-90 to 90 degrees}
%    lon = [1 x nobs] longitude {-180 to 360 degrees}
%    rtime = [1 x nobs] TAI time {double seconds since taiyear}
%    taiyear = [1 x 1] TAI year {1993 or 2000}
%
% Output:
%    efreq = [npts x 1] emissivity frequency points {wavenumber}
%    emis  = [npts x nobs] emissivity points {0 to 1, -999 if no data}
%

% Created: 21 Dec 2010, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlab/iasi/utils  % emis_iasi_DanZhou2
addpath /asl/matlab/science     % utc2tai & tai2utc

% The portion of Daniel Zhou (NASA Langley) database we have spans
% 2007/07 to 2008/06 with one file per month.  Interpolate across
% months ignoring year.
yyyymm = [200801, 200802, 200803, 200804, 200805, 200806, ...
          200707, 200708, 200709, 200710, 200711, 200712];
junk = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
center = (cumsum(junk) - junk/2)/sum(junk); % decimal year at month center
clear junk

% Emissivity frequency points to retain for output
efreq = [650,690,710,735,760,769,780,791,800,811,817,830,855,880,...
   905,930,945,960,985,1010,1025:10:1075,1090,1110,1135,1148,1159,1170,1180,...
   1205,1220,1230,1240,1260,1270,1280,1300:25:1350,1400:50:1650,...
   1690,1750,1775,1805,1835,1865,1895,1920:40:2000,2030,2040,2055, ...
   2070,2090,2105:10:2145,2160:10:2250,2270,2300,2380,2400,2430,...
   2460,2490,2515,2540:20:2580,2595,2610:20:2730,2755];
efreq = efreq'; %' [npts x 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inde = round(4*(efreq - 645) + 1);  % exact integers
nemis = length(inde);

if (nargin ~= 4)
   error('unexpected number of input arguments')
end
d = size(lat);
if (min(d) ~= 1 | length(d) ~= 2)
   error('unexpected dimensions for argument lat')
end
nobs = length(lat);
d = size(lon);
if (min(d) ~= 1 | length(d) ~= 2)
   error('unexpected dimensions for argument lon')
end
if (max(d) ~= nobs)
   error('arguments lat and lon must be the same length')
end
d = size(rtime);
if (min(d) ~= 1 | length(d) ~= 2)
   error('unexpected dimensions for argument rtime')
end
if (max(d) ~= nobs)
   error('argument rtime must be the same length as lat & lon')
end
if (min(rtime) <= 0)
   error('argument rtime contains values <= 0')
end
d = size(taiyear);
if (min(d) ~= 1 | max(d) ~= 1)
   error('unexpected dimensions for argument taiyear')
end
if (taiyear ~= 1993 & taiyear ~= 2000)
   error('taiyear must be 1993 or 2000')
end


% Convert rtime into decimal year
taistr = int2str(taiyear);
eval(['[minyear, month, day, dhour] = tai2utc' taistr '(min(rtime));'])
eval(['[maxyear, month, day, dhour] = tai2utc' taistr '(max(rtime));'])
maxyear = maxyear + 1;
nyears = maxyear - minyear + 1;
taistartofyear = zeros(1,nyears);
year = minyear;
for ii=1:nyears
   eval(['taistartofyear(ii) = utc2tai' taistr '(year,01,01,0);'])
   year = year + 1;
end
decimalyear = zeros(1,nobs);
for ii=1:(nyears-1)
   ip = find(rtime >= taistartofyear(ii) & rtime < taistartofyear(ii+1));
   decimalyear = (rtime(ip) - taistartofyear(ii)) / ...
      (taistartofyear(ii+1) - taistartofyear(ii));
end


% Declare output emis array
emis = zeros(nemis, nobs);


% Do early January (if any)
ip = find(decimalyear <= center(1));
np = length(ip);
if (np > 0)
disp('early jan')
   clo = center(12) - 1;
   chi = center(1);
   [elo] = emis_iasi_DanZhou2(yyyymm(12),lat(ip),lon(ip),inde);
   [ehi] = emis_iasi_DanZhou2(yyyymm(1 ),lat(ip),lon(ip),inde);
   ilo = find(elo(1,:) > 0);
   ihi = find(ehi(1,:) > 0);
   ionly = setdiff(ilo, ihi); % good elo & bad ehi
   ehi(:,ionly) = elo(:,ionly);
   ionly = setdiff(ihi, ilo); % good ehi & bad elo
   elo(:,ionly) = ehi(:,ionly);
   emis(:,ip) = ((ehi - elo)./(chi - clo)) .* ...
      (ones(nemis,1)*(decimalyear(ip) - clo)) + elo;
end


% Do late December (if any)
ip = find(decimalyear >= center(12));
np = length(ip);
if (np > 0)
disp('late dec')
   clo = center(12);
   chi = 1 + center(1);
   [elo] = emis_iasi_DanZhou2(yyyymm(12),lat(ip),lon(ip),inde);
   [ehi] = emis_iasi_DanZhou2(yyyymm(1 ),lat(ip),lon(ip),inde);
   ilo = find(elo(1,:) > 0);
   ihi = find(ehi(1,:) > 0);
   ionly = setdiff(ilo, ihi); % good elo & bad ehi
   ehi(:,ionly) = elo(:,ionly);
   ionly = setdiff(ihi, ilo); % good ehi & bad elo
   elo(:,ionly) = ehi(:,ionly);
   emis(:,ip) = ((ehi - elo)./(chi - clo)) .* ...
      (ones(nemis,1)*(decimalyear(ip) - clo)) + elo;
end


% Loop over other months
for ii=1:11
   clo = center(ii);
   chi = center(ii+1);
   ip = find(decimalyear >= clo & decimalyear <= chi);
   np = length(ip);
   if (np > 0)
%disp(['between months ' int2str(ii) ' and ' int2str(ii+1)])
      [elo] = emis_iasi_DanZhou2(yyyymm(ii  ),lat(ip),lon(ip),inde);
      [ehi] = emis_iasi_DanZhou2(yyyymm(ii+1),lat(ip),lon(ip),inde);
      ilo = find(elo(1,:) > 0);
      ihi = find(ehi(1,:) > 0);
      ionly = setdiff(ilo, ihi); % good elo & bad ehi
      ehi(:,ionly) = elo(:,ionly);
      ionly = setdiff(ihi, ilo); % good ehi & bad elo
      elo(:,ionly) = ehi(:,ionly);
      emis(:,ip) = ((ehi - elo)./(chi - clo)) .* ...
         (ones(nemis,1)*(decimalyear(ip) - clo)) + elo;
   end
end


%%% end of function %%%
