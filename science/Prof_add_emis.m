function [prof_out Qflag models ]=Prof_add_emis(prof_in, year, month, day, interp, mode, kind, only )
% function [prof_out Qflag models ]=Prof_add_emis(prof_in, year, month, day, interp, mode, kind, only )
%
% Add the Land emissivity to the corresponding land FoVs.
%
% prof_in : RTP profile containing the FoVs.
% year    : Calendar year of observation
% month   :          month
% day     :          day
%
%  ****** optinal arguments ******
%  (but if you use one, you must fill them all)
%
% interp  : Temporal interpolation (0, 1, 2)
%           0 - nearest month. (default)
%           1 - linear interpolation of 2 surrounding months
%           (not impl.) 2 - quadratic interpolation of 3 surrounding months
% mode    : Spacital smoothing (averaging)
%           'nearest' - get the nearest data pixel (default)
%           'averaged' - average inside of an AIRS FoVs (not coded yet)
% kind    : Kind of water emissivity function.
%            1 - cal_seaemis(satzen)   (default)
%                    - used for the uniform_clear
%            2 - cal_seaemis2(satzen, wspeed) 
%                    - same as above except includes wind speed
%            3 - cal_seaemiswu(satzen, wspeed) 
%                    - this is based on what JPL/AIRS uses.
%            4 - cal_Tdep_seaemis(satzem, wspeed, Tskin) 
%                    - experimental; no not use.
% only     : 'land' - change only land.
%          : 'sea'  - change only sea.
%          : 'all'  - chage all. (default)
%          : 'seaall' - change all as if it were Sea.
%
% OUTPUT:
% prof_out : is the same as prof_in, but with the adjusted emissivities.
% Qflag    : a quality indicator.
%
% Qflag = 0000  : Ok
%         0001  : Bad Land emiss
%         0010  : Over Land
%         0100  : Over Water  
%
% check for: Qflag==3
%
% models - string with the name of the emissivity models used.
%
%----------------------------------------------------------------------------
% Coastal areas:
%        The Wisconsin data (and the MODIS products used for it) 
% have no coastal data.  But, regardless, there are FoVs with fractional
% landfrac. We deal with this by adding the correcponding proportions 
% of land and sea emissivity.   
%
%----------------------------------------------------------------------------
% Ice:
% 	The Wisconsin data has no provision for frozen ocean.
% According to http://www.csgnetwork.com/h2ofreezecalc.html,
% the ocean water freezing point is 271.2K 
% (given salinity=35PSU and pressure of 1000 mbar) 
% We will test for lat>60N and stemp<271.2, and replace the ocean 
% emissivity by a "standard" Ice emissivity, taken from the South Pole.
%
%   B.I. - 2013.09.24
%
%----------------------------------------------------------------------------
% IRemis data files located at: /carrot/s1/imbiriba/SPECIALRTPFILES/iremis
%
% Breno Imbiriba - 2007.09.11

% Updated: 8 Aug 2011  Paul Schou - updated the dates to be 0 padded
%         18 Aug 2011  Paul Schou - added wrapTo180 around the rlon values
%          1 Sep 2011  Paul Schou - added nrho and rho calculations at the end of the file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 0. Setup

if(nargin==4)
  interp=0;
  mode='nearest';
  kind=1;
  only='all'
elseif(nargin~=8)
  error('You must provide 4 or 8 input arguments');
end 

% AIRS FOV diameter
dia=13.5;

%iremisdir='/carrot/s1/imbiriba/SPECIALRTPFILES/iremis';
%iremisdir=[ Sys_HOME() '/asl/iremis' ];
iremisdir='/asl/data/iremis/CIMSS/';
models='Land(Wisconsin:/asl/data/iremis/CIMSS/); ';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Find data file names:
  [iremis_files_ndate]=iremis_filedates();

  if(datenum(year,month,day)>=max(iremis_files_ndate) + 31)
    fprintf(['You requested year %d, but data base only goes till ' datestr(max(iremis_files_ndate)) '. Setting year to ' datestr(max(iremis_files_ndate),'yyyy') '.\n'],year);
    [year x x x x x]=datevec(max(iremis_files_ndate));
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1.1 Compute the matlab integer date number
  ndate=datenum([ num2str(year) '.' num2str(month) '.' num2str(day) ],'yyyy.mm.dd');
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1.2 Load the IREMIS file table (which is an index based on this matlab date number

  iremis_file_n=1:length(iremis_files_ndate);

  iclosest_date=interp1(iremis_files_ndate, iremis_file_n, ndate, 'linear','extrap');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1.3 - We want to be able to perform a quadratic fitting, if needed.
  %       We will save three dates, and three coeficients, for the interpolation:

  % Do I want two ahead or two behind?
  %idates=[nearest(iclosest_date)-1 nearest(iclosest_date) nearest(iclosest_date)+1]; 
  idates=[round(iclosest_date)-1 round(iclosest_date) round(iclosest_date)+1]; 

  % if I got out of range, set it to invalid -0-
  ioutofrange=find(idates<1 | idates>max(iremis_file_n));

  idates(ioutofrange)=0;
  if(all(idates==0))
    warning('Date is beyond limits. Using nearest.');
    interp=0;
    idates=[0, interp1(iremis_files_ndate, iremis_file_n, ndate, 'nearest','extrap'), 0];
  end   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1.4 Compute the mixing coeficients:
  if(interp>1)
    error('Sorry but only can do interp 0 or 1 for now');
  elseif(interp==1)
  %if(iclosest_date-nearest(iclosest_date)>=0)
  if(iclosest_date-round(iclosest_date)>=0)
    if(idates(3)~=0)
      % .   . * .
      %      x y
      x=1-(iclosest_date-idates(2));
      y=1-(-iclosest_date+idates(3));
      lin_coef=[0 x y];
    else
      % .   .   0
      %     x *
      x=1;
      lin_coef=[0 x 0];
      interp=0;
     end
  else
    if(idates(1)~=0)
      % . * .   .
      %  x y 
      x=1-(iclosest_date-idates(1));
      y=1-(-iclosest_date+idates(2));
      lin_coef=[x y 0];
    else
      % 0 * .   .
      %     x
      x=1;
      lin_coef=[0 x 0];
      interp=0;
    end
  end
  else
    if(idates(2)~=0)
      lin_coef=[0 1 0];
    elseif(idates(3)~=0)
      lin_coef=[0 0 1];
    else 
      lin_coef=[1 0 0];
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1.5 Load the required files:
  if(strcmp(only,'land') | strcmp(only,'all'))
    ineed=find(lin_coef~=0);
    for iload=1:numel(ineed);
      yyyy=datestr(iremis_files_ndate(idates(ineed(iload))),'yyyy');
      ddd=iremis_files_ndate(idates(ineed(iload)))-datenum(yyyy,'yyyy')+1;
      filedate=[yyyy num2str(ddd,'%03d')];
      filename=[iremisdir '/iremis_' filedate '.mat'];
      if(exist(filename,'file'))
	dat(iload)=load(filename);
      else
        disp(['Prof_add_emis: The emissivity file ' filename 'is missing! Will look for a file of another year...']);
        for itt1=[-1 1];
	  filedate=[num2str(str2num(yyyy)+itt1) num2str(ddd,'%03d')];
	  filename=[iremisdir '/iremis_' filedate '.mat'];
          if(~exist(filename))
            disp(['Attempt: file ' filename ' does not exist...']);
            filename='';
            continue;
          else
            disp(['Using file ' filename '.']);
            dat(iload)=load(filename);
            break  
          end
          % Find use an arbitrary existing file
          filedate='2006001';
	  filename=[iremisdir '/iremis_' filedate '.mat'];
          dat(iload)=load(filename);
        end
      end
      idxname=[iremisdir '/iremis_' filedate '_idx.mat'];
      idx(iload)=load(idxname);
    end
  end
  %ifovs=length(prof_in.rlat);

  loc(:,1)=prof_in.rlat; %(ifovs);
  loc(:,2)=wrapTo180(prof_in.rlon); %(ifovs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Get Emissivity:

wstr.landfrac=prof_in.landfrac;
switch kind
  case 1 
    wstr.satzen=prof_in.satzen;
    models=[models 'Water(cal_seaemis:' which('cal_seaemis') ')'];
  case 2
    wstr.satzen=prof_in.satzen;
    wstr.wspeed=prof_in.wspeed;
    models=[models 'Water(cal_seaemis2:' which('cal_seaemis2') ')'];
  case 3 
    wstr.satzen=prof_in.satzen;
    wstr.wspeed=prof_in.wspeed;
    models=[models 'Water(cal_seaemiswu:' which('cal_seaemiswu') ')'];
  otherwise 
  error('Prof_add_emis: not ready for kind=4');
end

if(strcmp(only,'land') | strcmp(only,'all'))
  if(interp==0)
     [nemis efreq emis Qflag] = Emis_get_all(loc, dat(1).land, idx(1).index_table, mode, kind, wstr, dia);
  elseif(interp==1)
     [nemis efreq emis Qflag] = Emis_get_all(loc, dat(1).land, idx(1).index_table, mode, kind, wstr, dia);
     [nemis2 efreq2 emis2 Qflag2] = Emis_get_all(loc, dat(2).land, idx(2).index_table, mode, kind, wstr, dia);
     emis = emis.*x + emis2.*y;
  end
else
  [nemis efreq emis]=Emis_get_water(loc, 'nearest', kind, wstr);
  Qflag=4*ones(size(loc(:,1)));
end

  % 2.1 Copy only the desired portion of data:
  if(strcmp(only,'land'))
    ifovs=find(wstr.landfrac>.99);
  elseif(strcmp(only,'sea'))
    ifovs=find(wstr.landfrac<.01);
  elseif(strcmp(only,'all'))
    ifovs=[1:length(wstr.landfrac)];
  elseif(strcmp(only,'seaall'))
    ifovs=[1:length(wstr.landfrac)];
  else
    warning(['Prof_add_emis: argument `only`==' only '. It must be land, sea, all, or seaall. Setting to `all`.']);
    ifovs=[1:length(estr.landfrac)];
  end



prof_out=prof_in;

prof_out.nemis(1,ifovs)=nemis(1,ifovs);

prof_out.emis=NaN(max(nemis), length(ifovs));
prof_out.efreq=NaN(max(nemis), length(ifovs));

% Copy emissivity to the profile structure.
% To speed up processes, we look for equal nemis, and do it as a chunk.

nemis_s=unique(nemis);
for i=1:length(nemis_s)
  inemis=find(nemis(ifovs(:))==nemis_s(i));
  prof_out.efreq(1:nemis_s(i),ifovs(inemis))=efreq(1:nemis_s(i),ifovs(inemis));
  prof_out.emis(1:nemis_s(i),ifovs(inemis))=emis(1:nemis_s(i), ifovs(inemis));
end


%------------------------------------------
% Test for ICE conditions:
% Iceman data:
% http://www.icess.ucsb.edu/modis/EMIS/html/ice.html

i_npole_ice = find(prof_in.stemp<271.3 & prof_in.rlat>60 & prof_in.landfrac<.001);
Modis_iceman01 = ...
	[0.9414 0.9324 0.9299 0.9295 0.9315 0.9364 0.9441 0.9541 0.9644 0.9771 0.9885 0.9945 0.9959 0.9932 0.9895 0.9872 0.9857 0.9858 0.9842 0.9835 0.9831 0.9827 0.9811 0.9799 0.9800 0.9800 0.9796 0.9805 0.9802 0.9807 0.9817 0.9802 0.9802 0.9795 0.9775 0.9772 0.9755 0.9741 0.9767 0.9783 0.9778 0.9796 0.9800 0.9813 0.9823 0.9836 0.9835 0.9838 0.9826 0.9828 0.9834 0.9846 0.9843 0.9832 0.9851 0.9849 0.9822 0.9827 0.9821 0.9803 0.9809 0.9792 0.9792 0.9780 0.9743 0.9752 0.9777 0.9786 0.9808 0.9791 0.9812 0.9807 0.9826 0.9803 0.9846 0.9834 0.9847 0.9829 0.9795 0.9817 0.9814 0.9812 0.9796 0.9780 0.9804 0.9777 0.9766 0.9776 0.9737 0.9720 0.9781 0.9724 0.9692 0.9724 0.9702 0.9677 0.9655 0.9617 0.9578 0.9547];
Modis_iceman01_efreq = 1000*...
	[0.6870 0.7103 0.7337 0.7570 0.7804 0.8037 0.8271 0.8504 0.8737 0.8971 0.9204 0.9438 0.9671 0.9905 1.0138 1.0372 1.0605 1.0838 1.1072 1.1305 1.1539 1.1772 1.2006 1.2239 1.2472 1.2706 1.2939 1.3173 1.3406 1.3640 1.3873 1.4106 1.4340 1.4573 1.4807 1.5040 1.5274 1.5507 1.5741 1.5974 1.6207 1.6441 1.6674 1.6908 1.7141 1.7375 1.7608 1.7841 1.8075 1.8308 1.8542 1.8775 1.9009 1.9242 1.9475 1.9709 1.9942 2.0176 2.0409 2.0643 2.0876 2.1109 2.1343 2.1576 2.1810 2.2043 2.2277 2.2510 2.2744 2.2977 2.3210 2.3444 2.3677 2.3911 2.4144 2.4378 2.4611 2.4844 2.5078 2.5311 2.5545 2.5778 2.6012 2.6245 2.6478 2.6712 2.6945 2.7179 2.7412 2.7646 2.7879 2.8113 2.8346 2.8579 2.8813 2.9046 2.9280 2.9513 2.9747 2.9980];

for ifov=1:numel(i_npole_ice)
  prof_out.nemis(:,i_npole_ice(ifov)) = 100;
  prof_out.efreq(1:100,i_npole_ice(ifov)) = Modis_iceman01_efreq;
  prof_out.emis(1:100,i_npole_ice(ifov)) = Modis_iceman01;
end
if(numel(i_npole_ice)>0)
  models = [models ' ice(Ts<271.3,Lat>60,Sea = Modis_iceman01)'];
end
%------------------------------------------


prof_out.nrho= prof_out.nemis;
prof_out.rho = (1.0 - prof_out.emis)/pi;


%for ic=1:length(ifovs)
%  prof_out.efreq(1:nemis(ifovs(ic)),ifovs(ic))=efreq(1:nemis(ifovs(ic)),ifovs(ic));
%  prof_out.emis(1:nemis(ifovs(ic)) ,ifovs(ic))= emis(1:nemis(ifovs(ic)),ifovs(ic));
%end

end




