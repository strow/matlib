cd ~/Desktop

start_time = datenum(2002,9,1);
stop_time  = datenum(2015,9,1);

d = read_netcdf_lls('co2_rta_aircraft-pfp_1_allvalid.nc');
co2 = d.value*1E6;
t = datenum(1970,1,1,0,0,double(d.time));

% Quick and dirty time for now
ndi = find( t >= start_time & t <= stop_time);
nd = length(ndi);
%if fit_type == 'rclr'
k = remove_6sigma(co2);
j = remove_6sigma(co2(k));
k = k(j);

t = t(ndi);
x = t-t(1);
y = co2(ndi);

[b stats] = Math_tsfit_lin_robust(x,y,4);
all_b = b;
all_rms = stats.s;
all_berr = stats.se;
all_bcorr = stats.coeffcorr;
%    all_resid(it(k),i) = stats.resid;
%    all_times(it(k),i) = fittime(k);
all_resid = stats.resid;
all_times = t;
% 
% 
% % Get lag-1 correlation (ignoring that we don't have all days)
% for i = 1:nf
%    y = squeeze(all_bt_resid(:,i));
%    k = remove_nan(y);
%    if length(k) > 100
%       l = xcorr(y(k),1,'coeff');
%       lag(i) = l(1);
%    else
%       lat(i) = NaN;
%    end
% end
% 
