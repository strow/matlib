function [] = fit_robust_one_lat_cldforcing(fin,fout,latid,fit_type,start_time,stop_time);

% Assumes  won't call with "robs" since already done
% Uses rclr for rcal, so bias = -(robs-rclr) is cloud forcing

addpath ~/Matlab/Math

latid_out = latid;

addpath /asl/matlib/aslutil
addpath ~/Matlab/Stats

load([fin int2str(latid)])
load_fairs

[nd,nf] = size(robs);
if nf ~= 2378
   disp('Error: robs wrong size');
   return;
end

all_b     = NaN(2378,10);
all_berr  = NaN(2378,10);
all_rms   = NaN(2378);
all_resid = NaN(nd,2378);
all_times = NaN(nd,2378);
all_bt_resid = NaN(nd,2378);
all_bcorr = NaN(2378,10,10);

% 3 basically means std(count over time) must be within 3x of a good channel
% keyboard
% 
btobs = real(rad2bt(f,robs'));
btcal = real(rad2bt(f,rclr'));
bias = btobs-btcal;
% Note minus so forcing is positive
bias = -bias';
% 

ig = goodchan_for_clear(count);

% In addition, remove 3-sigma (double pass) outliers for each channel
for i=ig
   k = remove_3sigma(robs(:,i));
   j = remove_3sigma(robs(k,i));
   ind(i).k = k(j);
end

% Quick and dirty time for now
dmtime = datenum(1958,1,1,0,0,rtime);

warning('off','all');
for i=ig
   it = ind(i).k;
   fittime = dmtime(it);

   switch fit_type
     case 'robs'
       fity = squeeze(robs(it,i));
     case 'rcal'
       fity = squeeze(rclr(it,i));
     case 'bias'
       fity = squeeze(bias(it,i));       
   end
   % Subset to specified times
   k = find( fittime >= start_time & fittime <= stop_time);
   x = fittime(k) - fittime(k(1));
   y = squeeze(fity(k));
   [b stats] = Math_tsfit_lin_robust(x,y,4);
   all_b(i,:) = b;
   all_rms(i) = stats.s;
   all_berr(i,:) = stats.se;
   all_bcorr(i,:,:) = stats.coeffcorr;
   all_resid(it(k),i) = stats.resid;
   all_times(it(k),i) = fittime(k);
end
warning('on','all');

% Get dr into dbt units
switch fit_type
  case {'robs', 'rcal'}
    deriv = drdbt(f,rad2bt(f,all_b(:,1)));
    dbt = all_b(:,2)./(deriv*1E3);
    dbt_err = all_berr(:,2)./(deriv*1E3);
  case 'bias'
    dbt = all_b(:,2);
    dbt_err = all_berr(:,2);
end

% Convert resid to BT units
switch fit_type
  case {'robs', 'rcal'}
    for i=1:nd
       all_bt_resid(i,:) = all_resid(i,:)./(deriv'*1E3);
    end
  case 'bias'
    for i=1:nd
       all_bt_resid(i,:) = all_resid(i,:);
    end
end

% Get lag-1 correlation (ignoring that we don't have all days)
for i = 1:2378
   y = squeeze(all_bt_resid(:,i));
   k = remove_nan(y);
   l = xcorr(y(k),1,'coeff');
   lag(i) = l(1);
end

clf;
plot(f(ig),dbt(ig));hold on;plot(f(ig),dbt_err(ig));title(['Lat: ' int2str(latid)]);
drawnow

% Convert resid to BT units
switch fit_type
  case {'robs', 'rcal'}
    save([fout int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');
%% Use if fitting clear calcs for cloudy random data
%    save([fout 'clrcal_' int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');
  case 'bias'
    save([fout int2str(latid_out)],'dbt','dbt_err','all*', 'lag');
end



