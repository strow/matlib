function [emis] = emis_iasi_DanZhou2(yyyymm,rlat,rlon,ind,loc);

% function [emis] = emis_iasi_DanZhou2(yyyymm,rlat,rlon,ind,loc);
%
% Return emissivity for IASI as determined by Dan Zhou (NASA Langley).
% Based on Dans "emis_sample_reader.f" given to L.Strow circa 2 August 2010
%
% Input:
%    yyyymm - [1 x 1] integer year*100 + month {200707 to 200806}
%    rlat   - [1 x nobs] latitude {degrees -90 to 90}
%    rlon   - [1 x nobs] longitude {degrees -180 to 360}
%    ind    - [npts x 1] desired IASI channel indices (in 1:8461)
%    loc    - Location of data file (optional dir string)
%             (default='/asl/data/iremis/DanZhou/Modified/')
%
% Output:
%    emis   - [npts x nobs] emissivity; -999 if no data
%

% Create: 11 Aug 2010, Scott Hannon - based on "emis_sample_reader.f" Dan Zhou
% Update: 14 Dec 2010, S.Hannon - update a few comments, no code changes
% Update: 21 Dec 2010, S.Hannon - version2 with "ind" argument created
% Update: 22 Dec 2010, S.Hannon - add missing "fclose" for ret file
% Update: 23 Dec 2010, S.Hannon - switch to "Modified" database
% Update: 24 Mar 2011, P.Schou - updated directory to /asl/data/iremis/DanZhou/
% Updade: 02 Jul 2013, B.Imbiriba - updatede location 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

nchan = 8461;    % number of IASI channels
nev = 10;        % number of eigenvectors
nrec = 720*360;  % number of records in ret file
% Note: grid points centered at lat=-89.75:0.5:89.75, lon=-179.75:0.5:179.75
nother = 5;      % number of other values preceeding coefs in ret file
% Note: "other" 1=landfrac; 2=lat; 3=lon; 4=surf pres; 5=surf T

% Data files
%datadir = '/home/strow/Transfer/DanZ/';
%datadir = '/asl/data/IASI/Emis_DanZhou/';
if(nargin()==5)
  datadir = loc;
else
  datadir = '/asl/data/iremis/DanZhou/Modified/';
end
evname = 'IASI_B_EV_FUNC_GLOBAL_V3.bin';  % eigenvector file name
retprefix = 'IASI_LAND_eFEOFA_';          % coef file prefix
retsuffix = '_0.5DEG_1SIGMA.bin';         % coef file suffix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input
if (nargin ~= 5)
   error('unexpected number of input arguments')
end
d = size(yyyymm);
if (length(d) ~= 2 | min(d) ~= 1 | max(d) ~= 1)
   error('unexpected dimesions for yyyymm')
end
if (yyyymm < 200707 | yyyymm > 200806)
   error('yyyymm is outside of allowed 200707-200806 range')
end
if (yyyymm > 200712 & yyyymm < 200801)
   error('yyyymm is outside of allowed 200707-200712,200801-200806 range')
end
d = size(rlat);
if (length(d) ~= 2 | min(d) ~= 1)
   error('unexpected dimensions for rlat')
end
nobs = max(d);
d = size(rlon);
if (length(d) ~= 2 | min(d) ~= 1)
   error('unexpected dimensions for rlon')
end
if (max(d) ~= nobs)
   error('rlat and rlon are not the same length')
end
d = size(ind);
if (length(d) ~= 2 | min(d) ~= 1)
   error('unexpected dimensions for ind')
end
if (max(d) > 8461)
   error('length of ind exceeds 8461')
end
if (min(ind) < 1 | max(ind) > 8461)
   error('ind contains values outside 1-8461 allowed range')
end
% Make sure ind is [npts x 1]
if (d(1) == 1)
   inde = ind'; %'
else
   inde = ind;
end
npts = length(inde);


% Determine filenames
inevf = [datadir evname];
ymstr = int2str(yyyymm);
rtvfile = [datadir retprefix ymstr retsuffix];


% Convert lat/lon to a record number
% Make sure (-90 <= lat < 90) and (-180 <= lon < 180)
xlat = rlat;
ii = find(xlat < -90);
xlat(ii) = -90;
ii = find(xlat >= 90);
xlat(ii) = 89.999; % lat must be less than 90
xlon = rlon;
ii = find(rlon >= 180); % lon must be less than 180
xlon(ii) = rlon(ii) - 360;
ii = find(xlon < -180);
xlon(ii) = -180;
%
n1 = floor( (180+xlon)/0.5 ) + 1;  % lon grid index 1-720
n2 = floor( ( 90+xlat)/0.5 ) + 1;  % lat grid index 1-360
nnn = (n1-1)*360 + n2; % index of ret file record

%%% this block for testing
%disp(['min(n1),max(n1)=',int2str(min(n1)) ', ' int2str(max(n1))])
%disp(['min(n2),max(n2)=',int2str(min(n2)) ', ' int2str(max(n2))])
%disp(['min(nnn),max(nnn)=',int2str(min(nnn)) ', ' int2str(max(nnn))])
%%%

% Read eigenvector file (no fortran record markers)
fid = fopen(inevf,'r','ieee-be');
emean = fread(fid,nchan,'real*8');
vcv = fread(fid,[nchan,nev],'real*8');
fclose(fid);


% Subset emean and vcv for ind
emean = emean(inde);
vcv = vcv(inde,:);


% Read emis retrieval coeffcient file
fid = fopen(rtvfile,'r','ieee-be');
n1 = nother + 1;
n2 = nother + nev;
ret = fread(fid,[n2,nrec],'real*4');
fclose(fid);


%%% this block for testing
%disp(['lat min&max:' num2str(min(ret(2,:))) ', ' num2str(max(ret(2,:)))])
%disp(['lon min&max:' num2str(min(ret(3,:))) ', ' num2str(max(ret(3,:)))])
%disp('all unique lat')
%unique(ret(2,:))
%disp('all unique lon')
%unique(ret(3,:))
%%%


% Determine indices of output with no data
ii = find(ret(1,:) < 0.2 & ret(6,:) == single(-999));
jj = find(ret(5,:) < 50);
ibad = union(ii,jj); % no data grid points indices
istat = zeros(1,nrec);
istat(ibad) = 1;
ibad = find(istat(nnn) == 1); % no data output indices
iok = setdiff(1:nobs,ibad);
nok = length(iok);


% Declare output
emis = -999*ones(npts,nobs);


% Compute emis
if (nok > 0)
   emis(:,iok) = emean*ones(1,nok);
   for iev = 1:nev
      emis(:,iok) = emis(:,iok) + vcv(:,iev)*ret(nother+iev,nnn(iok));
   end
   % Note sure why Dan does the following adjustment
   emis(:,iok) = 0.995 - exp(-0.693147181 - exp(emis(:,iok)));
end
 
return

%%% end of function %%%
