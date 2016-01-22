function [year, month, day, dhour] = tai2utc1993(tai);

%
% function [year, month, day, dhour] = tai2utc1993(tai);
%
% Convert TAI time to approximate UTC date & decimal hour
%
% Input:
%    tai : (1 x n) TAI time (seconds since 00:00 1 Jan 1993)
%
% Output:
%    year  : (1 x n) 4 digit integer year
%    month : (1 x n) 1 or 2 digit integer month
%    day   : (1 x n) 1 or 2 digit integer day
%    dhour : (1 x n) decimal hour
%

% Created: 9 July 2003 Scott Hannon
% Update: 07 July 2009, S.Hannon - fix rare bad day of month (eg 220838400
%    returned 0z 32 Dec 1999 instead of 0z 01 Jan 2000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Assign basic info

% Number of days per month
%        1  2  3  4  5  6  7  8  9 10 11 12
ndpm  =[31 28 31 30 31 30 31 31 30 31 30 31];  % non-leap year
ndpmly=[31 29 31 30 31 30 31 31 30 31 30 31];  % leap year

% Allowed years
allyears=1993:2099;
nall=length(allyears);

% Leap years every 4 years except century years unless divisible by 400
leapyears=1996:4:2099;


%%%
% Check input
if ( ndims(tai) > 2 | min(size(tai)) ~= 1 )
   error('tai must be a scaler or 1-D vector')
end
if (min(tai) < 0 | max(tai) > 3.2E+9)
   error('tai must be between 0 and 3.2E+9')
end
[irow,icol]=size(tai);
n=length(tai);



%%%
% Create TAI lookup tables for years and months
%
% Number of days in all allowed years
ndiy=365*ones(1,nall);
ii=find( ismember(allyears,leapyears) == 1);
ndiy(ii)=366;
%
% Seconds per day
spd=round(24*60*60); % exact integer
%
% Seconds since 1 Jan 1993 at start of each year
junk=[0, ndiy(1:(nall-1))]; % elapsed days since 1 Jan 1993
taiyear=round(cumsum(junk)*spd); % exact integer
%
% Seconds since start of year at start of each month
junk=[0, ndpm(1:11)]; % elapsed days of year at start of each month
taimonth=round(cumsum(junk)*spd); % exact integer
junk=[0, ndpmly(1:11)];
taimonthly=round(cumsum(junk)*spd); % exact integer


%%%
% Determine year, month, day, dhour
year=zeros(1,n);
month=zeros(1,n);
year=zeros(1,n);
year=zeros(1,n);
for ii=1:n
   % year
   ity=max( find(tai(ii) >= taiyear) );
   year(ii)=allyears(ity);
   junk=tai(ii) - taiyear(ity);
   ly=ismember(year(ii),leapyears);
   %
   % month
   if (ly == 1)
      month(ii)=max( find(junk >= taimonthly) );
      junk=junk - taimonthly( month(ii) );
   else
      month(ii)=max( find(junk >= taimonth) );
      junk=junk - taimonth( month(ii) );
   end
   %
   % day & dhour
   ned=floor( junk/spd ); % number of elapsed days in month
   junk=junk - spd*ned;
   day(ii)=1 + ned; % day of month
   dhour(ii)=junk/3600;
end


%%% end of function %%%
