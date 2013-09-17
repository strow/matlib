function [cpsize] = fake_cpsize(temp, iceflag, randomCpsize);

% function [cpsize] = fake_cpsize(temp, iceflag, randomCpsize);
%
% Generate a fake cloud particle size based on temperature.
% Values are randomized about the mean if randomCpsize=1.  Set
% iceflag=1 for ice and 0 for water.
%
% Input:
%    temp    = [1 x n] temperature {K}
%    iceflag = [1 x n] ice flag {1=true, 0=false(liquid water)}
%    randomCpsize= [1 x 1] optional randomization switch 
%        {1=on, 0=off(default), 20=fix water at 20 um, -1 = modis water dme climatology}
%
% Output:
%    cpsize  = [1 x n] particle size(diameter) {um}
%
%
% randomCpsize   = 1  ==> random dme_water (centered about 20 um), dme_ice = based on Scott Tcld, randomized
%                 20  ==> dme_water FIXED at 20 um, dme_ice = based on KN Liou Tcld, randomized
%                 -1  ==> dme_water based on MODIS climatology, dme_ice = based on Scott Tcld, randomized
%               +9999 ==> dme_water based on MODIS climatology, dme_ice = based on KN Liou Tcld, randomized

% Created: 05 Mar 2009, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ice size vs temp lookup table
icesize = [ 30,  90, 170];
icetemp = [213, 238, 263];
icestd  = [  4,  15,  30];
nicesize = length(icesize);
minicesize = 20;
maxicesize = 200;

% Water size vs temp lookup table
watsize = [ 15,   18,  21];
wattemp = [213,  243, 273];
watstd  = [  1,  1.5,   2];
nwatsize = length(watsize);
minwatsize = 13;
maxwatsize = 26;

% if randomCpsize == 0
%   watsize = [ 20,   20,  20];  %% keep water size at 20 um at all T
%   wattemp = [213,  243, 273];
%   watstd  = [  1,  1.5,   2];
%   nwatsize = length(watsize);
%   minwatsize = 13;
%   maxwatsize = 26;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input
if (nargin < 2)
   error('Insufficient input arguments')
end
if (nargin > 3)
   error('Too many input arguments')
end
if (nargin == 2)
   randomCpsize = 0;
end
d = size(temp);
if (length(d) ~= 2 | min(d) ~= 1)
   error('temp must be a [1 x n] array')
end
n = max(d);
d = size(iceflag);
if (length(d) ~= 2 | min(d) ~= 1)
   error('iceflag must be a [1 x n] array')
end
if (max(d) ~= n)
   error('iceflag must be the same length as temp')
end
if (nargin == 3)
   d = size(randomCpsize);
   junk = randomCpsize(1);
   if (length(d) ~= 2 | max(d) ~= 1 | length(intersect(junk,[0,1,20,-1])) ~= 1)
      error('randomCpsize must be a [1 x 1] scaler with value 0 or 1 or 20 or -1')
   end
end

% Determine indices of ice and water in iceflag
inda = 1:n;
indi = find(iceflag == 1);
nindi = length(indi);
indw = setdiff(inda,indi);
nindw = length(indw);

% Declare output
cpsize = zeros(1,n);

% Declare worjk arrays
cpstd  = zeros(1,n);
cprand = zeros(1,n);

% Process ice data (if any) : this is default setting randomCpsize = 0
if (nindi > 0)
   ilo = indi( find(temp(indi) <= icetemp(1)) );
   cpsize(ilo) = icesize(1);
   cpstd(ilo) = icestd(1);
   ihi = indi( find(temp(indi) > icetemp(nicesize)) );
   cpsize(ihi) = icesize(nicesize);
   cpstd(ihi) = icestd(nicesize);
   for ii=1:(nicesize-1)
      jj = indi( find(temp(indi) > icetemp(ii) & ...
                      temp(indi) <= icetemp(ii+1)) );
      cpsize(jj) = icesize(ii) + (temp(jj) - icetemp(ii))*...
         (icesize(ii+1)-icesize(ii))/(icetemp(ii+1)-icetemp(ii));
      cpstd(jj) = icestd(ii) + (temp(jj) - icetemp(ii))*...
         (icestd(ii+1)-icestd(ii))/(icetemp(ii+1)-icetemp(ii));
   end
end

% Process water data (if any) : this is default setting randomCpsize = 0
if (nindw > 0)
   ilo = indw( find(temp(indw) <= wattemp(1)) );
   cpsize(ilo) = watsize(1);
   cpstd(ilo) = watstd(1);
   ihi = indw( find(temp(indw) > wattemp(nwatsize)) );
   cpsize(ihi) = watsize(nwatsize);
   cpstd(ihi) = watstd(nwatsize);
   for ii=1:(nwatsize-1)
      jj = indw( find(temp(indw) > wattemp(ii) & ...
                      temp(indw) <= wattemp(ii+1)) );
      cpsize(jj) = watsize(ii) + (temp(jj) - wattemp(ii))*...
	 (watsize(ii+1)-watsize(ii))/(wattemp(ii+1)-wattemp(ii));
      cpstd(jj) = watstd(ii) + (temp(jj) - wattemp(ii))*...
         (watstd(ii+1)-watstd(ii))/(wattemp(ii+1)-wattemp(ii));
   end
end

% Apply randomization (if desired)
if (randomCpsize == 1)
   randn('state',sum(100*clock));
   junk = randn(1,n) .* cpstd;
   cpsize = cpsize + junk;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if randomCpsize == 20 | randomCpsize == 9999
  disp('>>>>> WARNING : reset ice dme using PCRTM settings')

  cpsizeXYZ = cpsize;
  if nindw > 0 & randomCpsize == 20
    cpsize(indw) = 20;   %% fixed particle water size
  end

  if nindi > 0
    %% need to fix ice particle sizes, according to PCRTM
    % coefficents from S-C Ou, K-N. Liou, Atmospheric Research
    % 35(1995):127-138.
    % for computing ice cloud effective size
    c0 = 326.3;
    c1 = 12.42;
    c2 = 0.197;
    c3 = 0.0012;

    tcld = temp(indi) - 273.16;
    oo = find(tcld < -60);
      tcld(oo) = -60;  % set a minimum, xianglei set -50
    oo = find(tcld > -20);
      tcld(oo) = -20;  % set a maximum, xianglei set -25                   
    cldde_ice = c0 + c1 * tcld + c2 * tcld.^2 + c3 * tcld.^3 ;
    cpsize(indi) = cldde_ice;
  end
  % plot(1:length(cpsize),cpsizeXYZ,1:length(cpsize),cpsize,'r')
  % disp('ret to cont'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply min/max check
ilo = indw(find(cpsize(indw) < minwatsize));
cpsize(ilo) = minwatsize;
ihi = indw(find(cpsize(indw) > maxwatsize));
cpsize(ihi) = maxwatsize;
ilo = indi(find(cpsize(indi) < minicesize));
cpsize(ilo) = minicesize;
ihi = indi(find(cpsize(indi) > maxicesize));
cpsize(ihi) = maxicesize;

%%% end of function %%%
