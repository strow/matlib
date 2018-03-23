function [] = fit_robust_one_lat(fin,fout,latid,fit_type,start_time,stop_time,inst);

smallsave = true;

addpath ~/Matlab/Math
addpath /asl/matlib/aslutil

latid_out = latid;

count = 1;  % Needed to diable new Matlab Mfile count (count string occurances)
            % Now can load in count in command below and not have problems later
load([fin int2str(latid)])

% Quick and dirty time for now
dmtime = datenum(1958,1,1,0,0,rtime);
dmtime = nanmean(dmtime,2);
ndi = find( dmtime >= start_time & dmtime <= stop_time);
nd = length(ndi);

if inst == 'cris'
   nf = 1305;
   load_fcris
   f = fcris;
   load /asl/matlib/cris/ch_std_from1317  % get ch_std_i
   count = count(ndi,1);
   dmtime = dmtime(ndi);
   robs = squeeze(nanmean(robs(ndi,ch_std_i,:),3));
   rcal = squeeze(nanmean(rcal(ndi,ch_std_i,:),3));
   btobs = real(rad2bt(f,robs'));
   btcal = real(rad2bt(f,rcal'));
   bias = btobs-btcal;
   bias = bias';
else
   % add robs, rcal, etc. for AIRS and IASI here.
   dmtime = dmtime(ndi);
   count = count(ndi,:);
end

if inst == 'iasi'
   nf = 8461;
   load_fiasi;
   f = fiasi;
   robs = robs(ndi,:);
   rcal = rcal(ndi,:);
   btobs = real(rad2bt(f,robs'));
   btcal = real(rad2bt(f,rcal'));
   bias = btobs-btcal;
   bias = bias';
end
   
if inst == 'airs'
   nf = 2378;
   load_fairs;
   robs = robs(ndi,:);
   rcal = rcal(ndi,:);
   btobs = real(rad2bt(f,robs'));
   btcal = real(rad2bt(f,rcal'));
   bias = btobs-btcal;
   bias = bias';
end

if inst == 'al1c'
   nf = 2645
   load_fairs;
%   load /home/sbuczko1/git/rtp_prod2/airs/util/sarta_chans_for_l1c.mat
%   keyboard
   f = fairs;
   robs = robs(ndi,:);
   rcal = rcal(ndi,:);
   btobs = real(rad2bt(f,robs'));
   btcal = real(rad2bt(f,rcal'));
   bias = btobs-btcal;
   bias = bias';
end

% 
% 
% if nf == 2378
%    load_fairs
% elseif nf == 8461
%    load_fiasi;
%    f = fiasi;
% elseif nf == 1305
% end

all_b        = NaN(nf,10);
pall_berr    = NaN(nf,10);
all_rms      = NaN(nf,10);
all_resid    = NaN(nd,nf);
all_bt_resid = NaN(nd,nf);
all_times    = NaN(nd,nf);
all_bcorr    = NaN(nf,10,10);

if nf == 2378
   ig = goodchan_for_clear(count);
elseif nf == 8461
   ig = 1:8461;
elseif nf == 1305
   ig = 1:1305;
elseif nf == 2645
   ig = 1:2645;
end

for i=ig
   k = remove_6sigma(robs(:,i));
   j = remove_6sigma(robs(k,i));
   ind(i).k = k(j);
end

warning('off','all');

for i=ig
   it = ind(i).k;
%   it = length(dmtime);it = 1:it;
   fittime = dmtime(it);

   switch fit_type
     case 'robs'
       fity = squeeze(robs(it,i));
     case 'rcal'
       fity = squeeze(rcal(it,i));
     case 'rclr'
       fity = squeeze(rclr(it,i));
     case 'bias'
       fity = squeeze(bias(it,i));       
   end
   x = fittime - fittime(1);
   y = squeeze(fity);
   [b stats] = Math_tsfit_lin_robust(x,y,4);
   all_b(i,:) = b;
   all_rms(i) = stats.s;
   all_berr(i,:) = stats.se;
   all_bcorr(i,:,:) = stats.coeffcorr;
   all_resid(it,i) = stats.resid;
   all_times(it,i) = fittime;
end
warning('on','all');

% Get dr into dbt units
switch fit_type
  case {'robs', 'rcal', 'rclr'}
    deriv   = drdbt(f,rad2bt(f,all_b(:,1)));
    dbt     = all_b(:,2)./(deriv*1E3);
    dbt_err = all_berr(:,2)./(deriv*1E3);
  case 'bias'
    dbt     = all_b(:,2);
    dbt_err = all_berr(:,2);
end

% Convert resid to BT units
switch fit_type
  case {'robs', 'rcal', 'rclr'}
    for i=1:nd
       all_bt_resid(i,:) = all_resid(i,:)./(deriv'*1E3);
    end
  case 'bias'
    for i=1:nd
       all_bt_resid(i,:) = all_resid(i,:);
    end
end

% Get lag-1 correlation (ignoring that we don't have all days)
for i = 1:nf
   y = squeeze(all_bt_resid(:,i));
   k = remove_nan(y);
   if length(k) > 100
      l = xcorr(y(k),1,'coeff');
      lag(i) = l(1);
   else
      lat(i) = NaN;
   end
end

if ~smallsave
   switch fit_type
     case {'robs', 'rcal', 'rclr'}
       save([fout int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');
       %% Use if fitting clear calcs for cloudy random data
%    save([fout 'clrcal_' int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');
     case 'bias'
       save([fout int2str(latid_out)],'dbt','dbt_err','all*','lag');
   end
else % smallsave option
% less output
   switch fit_type
     case {'robs', 'rcal','rclr'}
       save([fout int2str(latid_out)],'dbt','dbt_err','all_b','all_berr','all_bcorr','all_rms','lag');
     case 'bias'
       save([fout int2str(latid_out)],'dbt','dbt_err','all_b','all_berr','all_bcorr','all_rms','lag');
   end
end

% OLD saves
%    save([fout 'clrcal_' int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');


