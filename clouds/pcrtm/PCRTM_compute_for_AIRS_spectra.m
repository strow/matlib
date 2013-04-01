function [rad_allsky rad_clrsky tmpjunk] = PCRTM_compute_for_AIRS_spectra(nboxes,nlev, ncol, overlap, P, WCT, ICT, cc, TT, q, o3, Ps, Ts, sfctype,efreq,emis,zen_ang,co2, parname,ppath)

% this code is for simulating AIRS observation using changed PCRTM_V2.1 and
% using whatever given profiles as long as the pressure levels are fixed
% for each profile
% the viewing zenith angle can be from 0 to 70 degree

% nboxes         the number of profiles                                                                
% nlev           the number of layers, all boxes has same levels
% ncol           the number of subcolumn; if equal to -1, this is ONE COLUMN, CLOUD FRAC = 1
% overlap =1     maximum overlap
% overlap =2     random overlap
% overlap =3     maximum-random overlap
% WCT            water content in kg/kg
% ICT            ice content in kg/kg
% cc             cloud fraction profile, 0-1
% P              pressure profile in hPa
% TT             temperature profile in k
% q              specific humidity in g/kg
% o3             ozone mass mixing ratio in g/kg
% Ps             surface pressure in hPa
% Ts             surface temperature in k
% sfctype        the surface type of each boxes
% emis           the surface emissiviyt
% if sfctype<=0  emis is used, otherwise use default surface emissivity of 18 IGBP types
% zen_ang        viewing zenith angle in degree
% co2            co2 volume mixing ratio (ppmv)            
% parname        name of the input file for PCRTM
% ppath          the path where the PCRTM excutable file is

% error checking
if nargin ~= 20
  disp('wrong number of inputs');
  return;
end

%% sergio : originally all that was done is the file was opened for "append" so if it already existed, watch out!
eee = exist(parname);
if eee > 0
  rmer = ['!/bin/rm ' parname];
  eval(rmer)
end

% PCRTM fixed level pressure in mb (total 101 levels)
% T,q, o3 profiles need to be interpolated into these fixed levels
% cloud profiles not need to be interpolated, but need level_id which shows
% cloud are between Pres(level_id) and Pres(level_id +1)
Pres = [     0.0050,    0.0161,    0.0384,    0.0769,    0.1370,    ...
             0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    ...
             1.2972,    1.6872,    2.1526,    2.7009,    3.3398,    ...
             4.0770,    4.9204,    5.8776,    6.9567,    8.1655,    ...
             9.5119,   11.0038,   12.6492,   14.4559,   16.4318,    ...
            18.5847,   20.9224,   23.4526,   26.1829,   29.1210,   ...
            32.2744,   35.6505,   39.2566,   43.1001,   47.1882,   ...
            51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   ...
            77.2396,   83.2310,   89.5204,   96.1138,  103.0172,   ...
           110.2366,  117.7775,  125.6456,  133.8462,  142.3848,  ...
           151.2664,  160.4959,  170.0784,  180.0183,  190.3203,  ...
           200.9887,  212.0277,  223.4415,  235.2338,  247.4085,  ...
           259.9691,  272.9191,  286.2617,  300.0000,  314.1369,  ...
           328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  ...
           407.4738,  424.4698,  441.8819,  459.7118,  477.9607,  ...
           496.6298,  515.7200,  535.2322,  555.1669,  575.5248,  ...
           596.3062,  617.5112,  639.1398,  661.1920,  683.6673,  ...
           706.5654,  729.8857,  753.6275,  777.7897,  802.3714,  ...
           827.3713,  852.7880,  878.6201,  904.8659,  931.5236,  ...
           958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170,  ...
         1100.0000 ]; 

R = 287.05;
g = 9.806;

% coefficents from S-C Ou, K-N. Liou, Atmospheric Research
% 35(1995):127-138.
% for computing ice cloud effective size
c0 = 326.3;
c1 = 12.42;
c2 = 0.197;
c3 = 0.0012;
  
homepath = pwd;
                               
%  cloud assignment to each sub column wthin a gridbox
[unique_col_frac ucol_num ucol_num_same subcol_frac] = ...
  get_subcolumn_frac_v2(nboxes, nlev, abs(ncol), cc', overlap);
if ncol > 0
  tmpjunk.unique_col_frac = unique_col_frac;
  tmpjunk.ucol_num        = ucol_num;
  tmpjunk.ucol_num_same   = ucol_num_same;
  tmpjunk.subcol_frac     = subcol_frac;
%else
%  error('oooh have not figured out the ONE COLUMN, CFRAC = 1???  TEST CASE YET!')
end

% convert ozone from g/kg to ppmv
o3 = o3 *1e3/47.998*28.96;

ATMno = 6;

newT = zeros(1,101);
newh2o = zeros(1,101);
newozone = zeros(1,101);

mod_default_profile = load('Modtran_standard_profiles.mat');

def_P = mod_default_profile.Pres(:, ATMno);
def_T = mod_default_profile.T_z(:, ATMno +1);
def_h2o = mod_default_profile.h2o(:, ATMno);
def_o3 = mod_default_profile.o3(:, ATMno);
% convert default h2o from ppmv to g/kg
def_h2o  = def_h2o *1e-3 *18/28.96;

endsign=0;

tstart = tic;
iMOD0 = 25;
for ibox =1:nboxes
  if mod(ibox,iMOD0) == 0
    tcurrent = toc(tstart);
    fprintf(1,'processing %6i out of % 6i in avg %8.6f secs \n',ibox,nboxes,tcurrent/iMOD0);
    tstart = tic;
  end

  % use default profiles for Pres smaller than the smallest pressure of each
  % box, as the value is identical for each box, we use P(1) as the smallest pressure
  indU = find (P(1,ibox) >= Pres);
  % parameters for default models in upper atmosphere
  for ilev = 1: indU(end)
    newT(ilev) = interp1(log(def_P), def_T, log(Pres(ilev)));                    
    newh2o(ilev) = exp(interp1(log(def_P), log(def_h2o), log(Pres(ilev))));         
    newozone(ilev) = exp(interp1(log(def_P), log(def_o3), log(Pres(ilev)))); 
  end

  endsign = 0;
     
  % get profiles for Pres larger than or equal to P(1)
  %  2012, Feb 26 found a bug: specific humidity is zero at lower
  %  atmosphere in reanalysis
  idw = find(q(:,ibox)<=0.0);
  q(idw,ibox) = 1.0e-6;

  % get T, q, o3 profiles for lower atmosphere from reanalysis top to surface
  for ilev =indU(end) + 1:length(Pres)
    % for Pres below the bottom pressure, T, q, and o3 are set to the value at bottom pressure
    if (Pres(ilev)>= P(end,ibox))
      Pp = P(end,ibox);
    elseif (Pres(ilev)>= P(1,ibox) & Pres(ilev)< P(end,ibox))
      Pp = Pres(ilev);
    end

    %% sergio : orig was interp1(log(P),TT(:,ibox),log(Pp))
    newT(ilev) = interp1(log(P(:,ibox)), TT(:,ibox), log(Pp),[],'extrap');               
    %% sergio : orig was exp(interp1(log(P),q(:,ibox),log(Pp)))
    newh2o(ilev) = exp(interp1(log(P(:,ibox)), log(q(:,ibox)), log(Pp),[],'extrap'));    
    %% sergio : orig was exp(interp1(log(P),o3(:,ibox),log(Pp)))
    newozone(ilev) = exp(interp1(log(P(:,ibox)), log(o3(:,ibox)), log(Pp),[],'extrap')); 
  end
    
  % set a minimum cloud cover (0.001) to decide clear or cloudy over the box
  idc = find(cc(:,ibox)>0.001);  
        
  % convert pressure level to altitude level
  Z = zeros(nlev,1);
         
  Px = P(:,ibox);
  for ilev = length(Px)-1:-1:1
    Ttmp = (TT(ilev,ibox) + TT(ilev+1,ibox))/2;
    scaleH = R*Ttmp/g/1000;
    dz = scaleH*log(Px(ilev+1)/Px(ilev));
    Z(ilev) = Z(ilev+1) + dz; 
  end
      
  sub_opt(1:nlev) = 0.0;
  % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  % First, set cloud parameter over each gridbox
    
  % set cloud pressure as the midlevel pressure of the highest level, from Klein S.A. and Jakob C., 1999
  cldpres_layer = exp(0.5*(log(P(1:end-1,ibox)) + log(P(2:end,ibox))));  
         
  cldpres_layer(nlev) = P(nlev,ibox);
  % replace the last cloud layer pressure with surface pressure
  % May 27, 2012
  idp = find(P(:,ibox)>=Ps(ibox));
  if ~isempty(idp)
    cldpres_layer(idp(1)-1) = Ps(ibox);
  end

  % set a default cloud effective size and cloud type on each layer
  cldde_layer(1:nlev) = 40.0;
  cldphase_layer(1:nlev) =1;
         
  % get the corresponding id of P(:,ibox) to Pres
  %  P(i,ibox)>Pres(l) & P(i,ibox)<Pres(l+1) 
  %  this is for assigning cloud to the nearest level of atmosphere ( the fixed 101 levels in PCRTM)
  clear idp;
  for i=1 : length(cldpres_layer)
    itmp = find(cldpres_layer(i) < Pres);
    idp(i) = itmp(1);
  end
         
  %% remember from above   idc = find(cc(:,ibox)>0.001);  
  sergio_ice_opt(1:nlev)   = 0;
  sergio_water_opt(1:nlev) = 0;

  for i = 1: length(idc)
            
    ilev = idc(i);
    %############################################################
    %compute cirrus cloud effective size from cloud temperature,
    % cite from S-C Ou, K-N. Liou, Atmospheric Research
    % 35(1995):127-138.
    % temperature of ice clouds is the atmosphere temperature
    % cirrus cloud effective size are fitted in the range from - 20 to - 60C 
    % here the range are reduced for most cirrsus temperature
    % is in this range
    tcld = TT(ilev,ibox) - 273.16;  % temperature of ice clouds
    if tcld<-50
      tcld = -50;  % set a minimum 
    end
    if tcld>-25
      tcld = -25;  % set a maximum
    end
                      
    cldde_ice(ilev) = c0 + c1 * tcld + c2 * tcld^2 + c3 * tcld^3 ;
    cldde_liq(ilev) = 20; % Diameter in PCRTM

    % cloud phase and cloud effective size  in  micron
    if cldpres_layer(ilev) <= 440
      cldphase_layer(ilev) = 2;  % cirrus cloud
      cldde_layer(ilev) = cldde_ice(ilev);
    else
      cldphase_layer(ilev) = 1;  % water cloud
      cldde_layer(ilev) = cldde_liq(ilev);
    end
              
    if ICT(ilev,ibox)>1e-10
      % compute ice cloud optical depth from Ebert and Curry (1992, J. Geophys. Res., 
      %    vol. 97, pp. 3831-3836.
      qi = ICT(ilev,ibox)/ TT(ilev,ibox) *P(ilev,ibox)*100/R *1e3;  %change ice water content from kg/kg to g/m^3  
      ice_opt =  (0.003448 + 2.431/cldde_ice(ilev))*qi/cc(ilev,ibox)*(Z(ilev)-Z(min(ilev+1,nlev))) *1e3; 
        % * 0.5 * (P(ilev)+P(max(ilev-1,1)))/g *1e4;
    else
      ice_opt =0;
      qi = 0;
    end
                
    qw = WCT(ilev,ibox)/ TT(ilev,ibox) *P(ilev,ibox)*100/R *1e3;  %change liquid water content from kg/kg to g/m^3
    % ECMWF technical report, Klein S. A., 1999                   
    water_opt = 3 * qw/cldde_liq(ilev)/cc(ilev,ibox)*(Z(ilev)-Z(min(ilev+1,nlev))) *1e3; 

    sub_opt(ilev) = ice_opt + water_opt;  %%% <<<<<<<<<<<<<<<<<<-- should we do this??? SSM 03/28/2013
         
    sergio_ice_opt(ilev)   = ice_opt;
    sergio_water_opt(ilev) = water_opt;

  end
  %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %  Now, set cloud parameters over each sub-column
    
  % prepare all sky profiles of each sub-column for PCRTM  
  % only input one profile if some profiles are the same
  % ucol_num(ibox) is the number of unique profiles    
     
  cldpres = zeros(ucol_num(ibox),nlev);
  cldopt  = zeros(ucol_num(ibox),nlev);
  cldde   = zeros(ucol_num(ibox),nlev);
  cldphase = zeros(ucol_num(ibox),nlev);
  cld_id = zeros(ucol_num(ibox),nlev);
  cldnum = zeros(1,ucol_num(ibox));
    
  usrstr = ['box',int2str(ibox)];
  for icol =1:ucol_num(ibox)  

    % find cloud levels over unique sub-column
    idx = find (unique_col_frac(ibox,icol,:) ==1);
    % If this sub-column has clouds
    cldnum(icol) = length(idx);  % number of cloud layers
    if ~isempty(idx)
      for ic = 1: cldnum(icol)
        % cloud top pressure in hPa                   
        cldpres(icol,ic) = cldpres_layer(idx(ic));
        cldphase(icol,ic) = cldphase_layer(idx(ic));   
        cldde(icol,ic) = cldde_layer(idx(ic));
        cld_id(icol,ic) = idp(idx(ic));
        % cloud optical depth
        cldopt(icol,ic) = sub_opt(idx(ic));       
      end % cloud layers
    end
  end % unique cloud profile

  pcrtm_cloud_stats    %% sets most of "tmpjunk" here

%  ibotmpjunk.
%  tmpjunk    
%%%% gets dumped out in tmp.in as 
%%%%    0 = start cloud info
%%%%    number of levels
%%%%      [level cldtop OD dme cldphase]  .... repeated for the "N" levels for this cloud
%%%%    repeat for other M clouds
%%%%    0 = end cloud info
%  xcldpres = zeros(20,37);
%  tmpjunk.cldpres(ibox,:,:)  = xcldpres;
%  tmpjunk.cldphase(ibox,:,:) = xcldphase;
%  tmpjunk.cldde(ibox,:,:)    = xcldde;
%  tmpjunk.cldopt(ibox,:,:)   = xcldopt;
%  tmpjunk.cld_id(ibox,:,:)   = xcld_id;

  if ibox == nboxes  
    endsign=1;     
  end   

  make_PCRTM2AIRS_input(Ps(ibox), Ts(ibox),Pres, newT, newh2o, newozone,...
                        ucol_num(ibox),cldnum,cld_id, cldpres, cldopt, cldde, cldphase, ...
                        sfctype(ibox),efreq(:,ibox),emis(:,ibox),zen_ang(ibox), co2(ibox),...
                        parname, usrstr, endsign);
    
end

disp('ENDED MAIN LOOP 1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NewVSOrig version

iNewVSOrig = -1;    %% do the orig version; calls PCRTM_V2.1.exe_orig which defaults to reading in      pcrtm.in
iNewVSOrig = +1;    %% do the new  version; calls PCRTM_V2.1.exe      which reads in arbitrary located  pcrtm.in

if iNewVSOrig == -1

  % do the orig version; calls PCRTM_V2.1.exe_orig which defaults to reading in      pcrtm.in

  % change the path to run PCRTM
  cd(ppath);

  % delete old pcrtm.in and make new pcrtm.in
  unix(['rm -rf ', 'pcrtm.in']);
  fid = fopen([homepath,'/pcrtm.in.tmp'], 'r');
  vid = fopen('pcrtm.in', 'w');

  lala = pwd;
  fprintf(1,'sergio line 1 : pcrtm.in in %s \n',lala);
  fprintf(1,'sergio line 2 : tmp      in %s \n',[homepath,'/pcrtm.in.tmp'])

  for iline = 1:5
    str = fgetl( fid);
    fprintf(vid, '%s\n', str);
    fprintf(1, '%s\n', str);
  end
  fprintf(vid, '%d\n', 1);

  fprintf(1, '%s\n', [native2unicode(39),parname,native2unicode(39)]);
  fprintf(vid, '%s\n', [native2unicode(39),parname,native2unicode(39)]);

  fclose(vid);
  fclose(fid);

  disp('CALLING PCRTM')
  system('./PCRTM_V2.1.exe_orig');

elseif iNewVSOrig == +1

  % do the new version; calls PCRTM_V2.1.exe_orig which defaults to reading in      pcrtm.in

  % change the path to run PCRTM
  cd(ppath);

  lalaNAME1 = [parname '.redirect_stdin'];
  lalaNAME2 = [parname '.pcrtm.in'];
  fprintf(1,'  name of redirect stdin file     is %s \n',lalaNAME1);
  fprintf(1,'  which will ask PCRTMexe to read in %s \n',lalaNAME2);

  % delete old pcrtm.in and make new pcrtm.in
  unix(['rm -rf ',lalaNAME1]);
  unix(['rm -rf ',lalaNAME2]);

  if ~exist([homepath,'/pcrtm.in.tmp'],'file')
    fprintf(1,'looking for %s \n',[homepath,'/pcrtm.in.tmp']);
    error('pcrtm.in.tmp file not found');
  end

  fid = fopen([homepath,'/pcrtm.in.tmp'], 'r');
  vid = fopen(lalaNAME2, 'w');

  lala = pwd;
  fprintf(1,'sergio line 1 : pcrtm.in in %s \n',lala);
  fprintf(1,'sergio line 2 : tmp      in %s \n',[homepath,'/pcrtm.in.tmp'])

  for iline = 1:5
    str = fgetl( fid);
    fprintf(vid, '%s\n', str);
    fprintf(1, '%s\n', str);
  end
  fprintf(vid, '%d\n', 1);

  fprintf(1, '%s\n', [native2unicode(39),parname,native2unicode(39)]);
  fprintf(vid, '%s\n', [native2unicode(39),parname,native2unicode(39)]);

  fclose(vid);
  fclose(fid);

  str = ['''' lalaNAME2 ''''];
  fidjunk = fopen(lalaNAME1,'w');
  fprintf(fidjunk,str);
  fclose(fidjunk)

  disp('CALLING PCRTM')
  runner = ['!./PCRTM_V2.1.exe < ' lalaNAME1];
  eval(runner)
end

if iNewVSOrig == +1
  rmer = ['!/bin/rm ' lalaNAME1 ' ' lalaNAME2];
  eval(rmer)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%Prof. Huang  revised%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
disp('RECONSTRUCTING RADS')
% coefficients database used for compute spectral Jacob 
% the PC loadings are precomputed with fixed size
fid=fopen(['../InputDir/Pc60_data_AIRS.dat'],'rb','l');

% data structure 
% num of bands, number of PCs in each band, and number of channels in each band
% these are all fixed value, reading them just for double-check
numbnd		= fread(fid, 1, 'float'); % should equal to 3
numPC 		= fread(fid, numbnd, 'float'); % should be [60, 60, 60]
					  % 60 PC in each band
numch		= fread(fid, numbnd, 'float'); %should be [1262, 602, 514]
					  % number of channels in each band 
					  % in total 2378 channels
						
% 1500 is used here for no numch >1500. Not all 1500 columns to be used
PCcoef = zeros(numbnd,numPC(1),1500)+NaN; % Each band has 60 PC loadings, 
				 % each loading has a dimension of numch (<1500)
Pstd = zeros(numbnd,1500)+NaN;

IDX = cell(numbnd, 1);

for i = 1:numbnd
  IDX{i} = 1:numch(i);
  if i >1
    IDX{i} = IDX{i} + sum(numch(1:i-1));
  end
end

for ib =1:numbnd
  for ip =1:numPC(ib)
    PCcoef(ib, ip, 1:numch(ib)) = reshape(fread(fid, numch(ib), 'float'), 1, 1, numch(ib));
  end
  Pstd(ib,1:numch(ib))= reshape(fread(fid, numch(ib), 'float'), 1, numch(ib));
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
cd (homepath);
     
if exist([parname,'.out'], 'file')
  tmp = load([parname,'.out']);
  PClen = sum(numPC);
    
  count = 0;
  PC_allsky =0;
  PC_clrsky =0;  
  for ibox = 1: nboxes
    % for each box, there are (PCend - PCstart)= PClen * (ucol_num(ibox) +1) PC scores  
    PCstart = count * PClen +1;
    count = count + ucol_num(ibox) +1;
    PCend = count * PClen ;

    %[ibox PClen PCstart PCend]
    %whos tmp ucol_num
                 
    tmprad = reshape(tmp(PCstart:PCend, 2), PClen,ucol_num(ibox) +1);
        
    % calculate average PCs of the box for all condition       
    rad=zeros(PClen,1);
    for j = 1:ucol_num(ibox)
      rad = rad +tmprad(:,j) .* ucol_num_same(ibox,j);
    end
    PC_allsky = reshape(rad/abs(ncol),numPC(1),numbnd);
        
    % clear sky spectra
    PC_clrsky = reshape(tmprad(:, ucol_num(ibox) +1),numPC(1), numbnd );  
        
    % convert PC domain to channel domain    
    
    for ib =1:numbnd  % band 1 to 3
      % for all-sky
      aa = PC_allsky(:, ib)' * reshape(PCcoef(ib, :, 1:numch(ib)), numPC(ib), numch(ib));		  
      rad_allsky(IDX{ib},ibox) = ((aa +1) .* Pstd(ib,1:numch(ib)))'*1e-3; % W per m^2 per sr per cm^-1

      %  for clear-sky
      bb = PC_clrsky(:, ib)' * reshape(PCcoef(ib,:,1:numch(ib)), numPC(ib), numch(ib)); 
      rad_clrsky(IDX{ib},ibox) = ((bb+1) .* Pstd(ib,1:numch(ib)))'*1e-3;  % W per m^2 per sr per cm^-1
    end  % end of bands
  end    % end of boxes
end      % end of output file


if exist('rad_allsky') & exist('rad_clrsky')
  plot(1:length(tmpjunk.totalODice),rad_allsky(1291,:),'bs-',...
       1:length(tmpjunk.totalODice),rad_clrsky(1291,:),'kx-')
  title('PCRTM cals for Rad 1231 : (b) cloudy (k) clear sky')
  pause(0.1)
end


