function [] = fit_robust_one_lat_cris(fin,fout,latid,fit_type,start_time,stop_time);

smallsave = false;

addpath ~/Matlab/Math
addpath /asl/matlib/aslutil

latid_out = latid;

load([fin int2str(latid)])

nf = 1305
load_fcris
f = fcris;

% Quick and dirty time for now
dmtime = datenum(1958,1,1,0,0,nanmean(rtime,2));
ndi = find( dmtime >= start_time & dmtime <= stop_time);
nd = length(ndi);

dmtime = dmtime(ndi);
robs = robs(ndi,:,:);
%rcal = rcal(ndi,:);
%rclr = rclr(ndi,:);
count = count(ndi,:);


all_b        = NaN(nf,10,9);
all_berr     = NaN(nf,10,9);
all_rms      = NaN(nf,10,9);
all_resid    = NaN(nd,nf,9);
all_bt_resid = NaN(nd,nf,9);
all_times    = NaN(nd,nf,9);
all_bcorr    = NaN(nf,10,10,9);

for ifov = 1:9
btobs(:,:,ifov) = real(rad2bt(f,squeeze(robs(:,:,ifov))'));
end
% btcal = real(rad2bt(f,rcal'));
% bias = btobs-btcal;
% bias = bias';

if nf == 2378
   ig = goodchan_for_clear(count);
elseif nf == 8461
   ig = 1:8461;
elseif nf == 1305
   ig = 1:1305;
end

% In addition, remove 3-sigma (double pass) outliers for each channel
% for i=ig
%    k = remove_6sigma(robs(:,i));
%    j = remove_6sigma(robs(k,i));
%    ind(i).k = k(j);
% end

%if fit_type == 'rclr'
for ifov = 1:9

for i=ig
      k = remove_6sigma(robs(:,i));
      j = remove_6sigma(robs(k,i));
      ind(i).k = k(j);
   end
   %end

warning('off','all');
for i=ig
%for i=1569
it = ind(i).k;
   fittime = dmtime(it);

   switch fit_type
     case 'robs'
       fity = squeeze(robs(it,i,ifov));
     case 'rcal'
       fity = squeeze(rcal(it,i));
     case 'rclr'
        fity = squeeze(rclr(it,i));
     case 'bias'
       fity = squeeze(bias(it,i));       
   end
   % Subset to specified times
%   k = find( fittime >= start_time & fittime <= stop_time);
%   x = fittime(k) - fittime(k(1));
   x = fittime - fittime(1);
%   y = squeeze(fity(k));
   y = squeeze(fity);
   [b stats] = Math_tsfit_lin_robust(x,y,4);
   all_b(i,:,ifov) = b;
   all_rms(i,ifov) = stats.s;
   all_berr(i,:,ifov) = stats.se;
   all_bcorr(i,:,:,ifov) = stats.coeffcorr;
%    all_resid(it(k),i) = stats.resid;
%    all_times(it(k),i) = fittime(k);
   all_resid(it,i,ifov) = stats.resid;
   all_times(it,i,ifov) = fittime;
end
warning('on','all');

% Get dr into dbt units
switch fit_type
  case {'robs', 'rcal', 'rclr'}
    deriv   = drdbt(f,rad2bt(f,all_b(:,1,ifov)));
    dbt     = all_b(:,2)./(deriv*1E3);
    dbt_err = all_berr(:,2,ifov)./(deriv*1E3);
  case 'bias'
    dbt     = all_b(:,2);
    dbt_err = all_berr(:,2);
end

% Convert resid to BT units
switch fit_type
  case {'robs', 'rcal', 'rclr'}
    for i=1:nd
       all_bt_resid(i,:,ifov) = all_resid(i,:,ifov)./(deriv'*1E3);
    end
  case 'bias'
    for i=1:nd
       all_bt_resid(i,:) = all_resid(i,:);
    end
end

% Get lag-1 correlation (ignoring that we don't have all days)
for i = 1:nf
   y = squeeze(all_bt_resid(:,i,ifov));
   k = remove_nan(y);
   if length(k) > 100
      l = xcorr(y(k),1,'coeff');
      lag(i,ifov) = l(1);
   else
      lat(i,ifov) = NaN;
   end
end

end


% clf;
% plot(f(ig),dbt(ig));hold on;plot(f(ig),dbt_err(ig));title(['Lat: ' int2str(latid)]);
% drawnow

% Convert resid to BT units

if ~smallsave
   switch fit_type
     case {'robs', 'rcal', 'rclr'}
       save([fout int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');
       %% Use if fitting clear calcs for cloudy random data
%    save([fout 'clrcal_' int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');
     case 'bias'
       save([fout int2str(latid_out)],'dbt','dbt_err','all*', 'lag');
   end
else % smallsave option
% less output
   switch fit_type
     case {'robs', 'rcal','rclr'}
       save([fout int2str(latid_out)],'dbt','dbt_err');
       %% Use if fitting clear calcs for cloudy random data
%    save([fout 'clrcal_' int2str(latid_out)],'dbt','dbt_err','all*', 'lag','deriv');
     case 'bias'
       save([fout int2str(latid_out)],'dbt','dbt_err');
   end
end



