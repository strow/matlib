function [unique_col_frac,ucol_num,ucol_num_same,subcol_frac] = get_subcolumn_frac_v2(nboxes, nlev, ncol, cc, overlap)

% this subroutinue is to assign cloud or clear to each sub-column of each of the "nboxes"
% box by threecloud overlap schemes
%   overlap = 1 maximum overlap
%           = 2 random overlap
%           = 3 maximum-random overlap
%   nboxes is the number of AIRS fovs being processed
%   nlev   is the number of layers of cloud profiles
%   ncol   is the number of sub-column
%   cc     is the cloud fraction profile (nboxes X nlev)
% so eg for 60 lev ERA, 100 AIRS fovs, and 50 subcolumns we have get_subcolumn_frac_v2(100, 60, 50, cc, 3)
%
% output
%   subcol_frac                    nboxes * ncol * nlev , the cloud assignments
%   unique_col_frac                nboxes * ucol_num * nlev, unique cloud assignments
%   ucol_num(nboxes)               number of unique cloud assignments
%   ucol_num_same(nboxes,ucol_num) the number of same cloud assignments for each unique cloud assignment
% with
%   sum(ucol_num_same(iboxes,i), i=1:ucol_num) = ncol   : plot(sum(ucol_num_same')) = ncol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% so basically ncol     = number of subgrid points we want to simulate for each of the ii grid points
%%              cc(:,ii) = ERA cloud fraction for the ii th grid point
%%              subcol_frac(ii,:,:) = subcol_frac(ii,ncol,nlevs) = random instances of putting clouds, 
%%                                    given info in cc(:,ii). BUT MANY ARE DUPLICATES OF EACH OTHER
%%              ucol_num(ii)         = tells you how many UNIQUE subcolumns there are, one of which is the ALL CLEAR case
%%              ucol_num_same(ii,:)  = tells how many unique subcolumns there are, for each of the above
%%              unique_col_frac(ii,:,:) = unique_col_frac(ii,X<50,nlevs) has 50 as second dimension, 
%%                                        but only fill upto X since the subcol_fracs are NOT all unique
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mm,nn] = size(cc);
if mm ~= nboxes
  error('mm ~= nboxes')
end

idd = find(cc<0.001);
cc(idd) = 0.0;

% add a clear layer over the top of atmosphere   
temp_cc(:,2:nlev+1) = cc;
temp_cc(:,1)        = 0.0;

boxpos = zeros(nboxes,ncol);   %% sergio
%% so this is basically delta : delta : 1 where delta = 1/ncol
for icol=1:ncol      
  boxpos(:,icol) = (icol - 0.5)/ncol;  %% so go between 0.5/ncol and (1 - 0.5/ncol)
end

% initialize the sub-column with 0 (clear sky)
subcol_frac = zeros(nboxes,ncol,nlev);
cloud_threshold = zeros(nboxes,ncol);
for ilev = 1:nlev
  if ilev == 1
    if overlap == 1
      % maximum overlap
      cloud_threshold = boxpos;
    else
      % random overlap or maximum and random overlap
      cloud_threshold = rand(nboxes,ncol);
    end
  end

  ran = rand(nboxes,ncol);
  for icol=1:ncol
    if overlap == 1
      % Maximum overlap
      tmin(1:nboxes,icol) = 0.0;
      flag(1:nboxes,icol) = 1;
    elseif overlap == 2
      % Random overlap
      tmin(1:nboxes,icol) = 0.0;
      flag(1:nboxes,icol) = 0;
    elseif overlap == 3
      % Maximum and Random overlap
      for ibox=1:nboxes
        tmin(ibox,icol) = min(temp_cc(ibox,ilev), temp_cc(ibox,ilev+1));
        %% if cloud_threshold(ibox,icol) < min(temp_cc(ibox,ilev), temp_cc(ibox,ilev+1)) & cloud_threshold(ibox,icol)>0
        if cloud_threshold(ibox,icol) > 0 & cloud_threshold(ibox,icol) < tmin(ibox,icol)
          flag(ibox,icol) = 1;
        else
          flag(ibox,icol) = 0;
        end
      end
    end
  
    %ranx=get_ran(nboxes, seed);
    %% CT(ilev+1) --> CT(ilev) * flag +  (1-flag)*(tmin + (1-tmin)*ran))    
    %% so for max,     cloud_threshold(ilev+1,icol) -> cloud_threshold(ilev,icol) == boxpos
    %%    for random,  cloud_threshold(ilve+1,icol) -> ran
    %%    for MR,      cloud_threshold(ilev+1,icol) -> cloud_threshold(ilev,icol) * flag +  (1-flag)*(tmin + (1-tmin)*ran))
    cloud_threshold(:,icol) = cloud_threshold(:,icol) .*flag(:,icol) + ...
                    (1-flag(:,icol)) .* (tmin(:,icol) + (1 - tmin(:,icol)) .* ran(:,icol));
  end
      
  % effectively code below is (ibox = 1)
  % for icol = 1:ncol
  %   if cc(ilev) > cloud_threshold(icol)
  %     subcol_frac(icol, ilev) = 1;
  %   else
  %     subcol_frac(icol, ilev) = 0;
  %  end
  %
  % (a) for MAXIMUM model this maps to following :
  %     suppose cc       looks as shown in left panel
  %     and     boxpos   looks as shown in right panel; note is has a minimum value of delta
  %
  %         ^                                   ^
  %       1 |                             1.0   |          /
  %  ilev 2 |                                   |         /
  %       . |                                   |        /
  %       . |                              cc2  |-------/
  %       M |---                                |      /|          Slope ~ 1/Ncol
  %         |   |                          cc1  |-----/ |          Which means to get to ccx you need to
  %       P |---------                          |    /| |            have gone ccx * Ncol boxes
  %         |        |                          |   / | |
  %       Q |---------                          |  /  | |
  %         |                             delta |-/   | |
  %       N +------------------+>               +-----+-+------+->
  %         0                 1.0               1   N1   N2      Ncol    
  %
  % then sub_colfrac1 will look exactly like plot on left!!!, except horizontal axis will be between 0 and Ncol
  % Reason :  for levs 1:M, cc(i) = 0 is always less than boxpos for all columns < Ncol (since min(boxpos) = delta).
  %           so those subcol_frac are 0
  %           levs M : Q have cc >= cc1 ==> for the N1 columns where cc1 > boxpos, these will get subcol_frac(icol,ilev) = 1
  %                      AND THIS WILL STRETCH from icol = 1 to icol = I1 where I1 = cc1 * Ncol!!!!!!!
  %           levs P : Q have cc = cc2 ==> for the N2 columns where cc2 > boxpos, these will get subcol_frac(icol,ilev) = 1
  %                      AND THIS WILL STRETCH from icol = 1 to icol = I2 where I2 = cc2 * Ncol!!!!!!!
  %           for levs Q:N, cc(i) = 0 is always less than boxpos for all columns < Ncol (since min(boxpos) = delta).
  %           so those subcol_frac are 0
  %
  %         ^                                   ^
  %       1 |                             1.0   |          /
  %  ilev 2 |                                   |         /
  %       . |                                   |        /
  %       . |                              cc2  |-------/
  %       M |---                                |      /|          Slope ~ 1/Ncol
  %         |   |                          cc1  |-----/ |          Which means to get to ccx you need to
  %       P |---------                          |    /| |            have gone ccx * Ncol boxes
  %         |        |                          |   / | |
  %       Q |---------                          |  /  | |
  %         |                             delta |-/   | |
  %       N +----+---+---------+>               +-----+-+------+->
  %         0   N1   N2       Ncol               1  N1   N2     Ncol  
  %
  % (b) for RANDOM at each level, RH plot is replaced   ^
  %     by random numbers (between 0 and 1)         1.0 |                                   +
  %                                                     |                     + +
  %                                                     |      +                 
  %                                                 cc2 |-------------------------------+---------------
  %                                                     |      
  %                                                     |   +
  %                                                 cc1 |-----+-----------------------------------------
  %                                                     |                 
  %                                                     | +               +
  %                                                     |              +
  %                                                 0.0 +----------------------------------------+--->
  %                                                      1                                      Ncol
  %
  % then sub_colfrac1 will look rather random, starting below level M
  % Reason :  for levs 1:M, cc(i) = 0 is almost always less than random number for all columns < Ncol
  %           so those subcol_frac are 0
  %           levs M : Q have cc >= cc1 ==> any of the subcols which have cloud_threshold < cc1 will get populated
  %           levs P : Q have cc =  cc2 ==> any of the subcols which have cloud_threshold < cc2 will get populated
  %           for levs Q:N, cc(i) = 0 is always less than boxpos for all columns < Ncol
  %           so those subcol_frac are 0
  %
  % (c) for MAX RAN, this is a beast!! but seems to work because of tmin = min(cc(i),cc(i+1)) :
  %        if there is no cloud at any of those two then tmin becomes 0,   and essentially flag --> 0 which means we get the RANDOM OVERLAP
  %        if there is    cloud at any of those two then tmin becomes > 0, and essentially flag --> 1 which means we get the MAXIMUM OVERLAP  
  %  
  % cloud assignment
  for icol = 1:ncol
    for ibox = 1:nboxes
      if cc(ibox,ilev) > cloud_threshold(ibox,icol)
        subcol_frac(ibox, icol, ilev) = 1;
      else
        subcol_frac(ibox, icol, ilev) = 0;
      end
    end
  end  
end

unique_col_frac = zeros(nboxes,ncol,nlev);
for ibox=1:nboxes
  subcol = reshape(subcol_frac(ibox,:,:),ncol,nlev);
  [fracnew,pos1,pos2] = unique(subcol, 'rows');
      
  ucol_num(ibox) = length(pos1); % the number of unique cloud profiles
      
  unique_col_frac(ibox,1:ucol_num(ibox),:) = fracnew;
      
  % the number of same cloud profiles for each unique cloud profile
  for i=1:ucol_num(ibox)
    ucol_num_same(ibox,i) = length(find (ismember(subcol, fracnew(i,:),'rows') ==1 ));
  end
end
