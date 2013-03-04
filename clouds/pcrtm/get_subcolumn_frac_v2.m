function [unique_col_frac,ucol_num,ucol_num_same,subcol_frac] = get_subcolumn_frac_v2(nboxes, nlev, ncol, cc, overlap)

% this subroutinue is to assign cloud or clear to each sub-column of each of the "nboxes"
% box by threecloud overlap schemes
%   overlap = 1 maximum overlap
%           = 2 random overlap
%           = 3 maximum-random overlap
%   nlev is the number of layers of cloud profiles
%   ncol is the number of sub-column
%   cc is the cloud fraction profile
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

idd = find(cc<0.001);
cc(idd) = 0.0;

% add a clear layer over the top of atmosphere   
tcc(:,2:nlev+1) = cc;
tcc(:,1) = 0.0;

boxpos = zeros(nboxes,ncol);   %% sergio
for icol=1:ncol      
  boxpos(:,icol) = (icol - 0.5)/ncol;
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
    if overlap ==1
      % Maximum overlap
      tmin(1:nboxes,icol) = 0.0;
      flag(1:nboxes,icol) = 1;
    elseif overlap ==2
      % Random overlap
      tmin(1:nboxes,icol) = 0.0;
      flag(1:nboxes,icol) = 0;
    elseif overlap ==3
      % Maximum and Random overlap
      for ibox=1:nboxes
        tmin(ibox,icol) = min(tcc(ibox,ilev), tcc(ibox,ilev+1));
        if cloud_threshold(ibox,icol)<min(tcc(ibox,ilev), tcc(ibox,ilev+1)) & cloud_threshold(ibox,icol)>0
          flag(ibox,icol) = 1;
        else
          flag(ibox,icol) = 0;
        end
      end
    end
  
    %ranx=get_ran(nboxes, seed);     
    cloud_threshold(:,icol) = ...
      cloud_threshold(:,icol) .*flag(:,icol) + ...
      (1-flag(:,icol)) .* (tmin(:,icol) + (1 - tmin(:,icol)) .* ran(:,icol));
  end
      
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
