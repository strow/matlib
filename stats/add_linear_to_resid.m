% Need BT

g = load('statlat30');
bt = real(rad2bt(f,g.robs'));

for j=1:2378
   k = find(~isnan(all_times(:,j)),1);
   if length(k) > 0
      ts(j) = k;
   else
      ts(j) = NaN;
   end
end

[nd,~] = size(all_resid);

for j=1:2378
   if ~isnan(ts(j))
      x = all_times(:,j) - all_times(ts(j),j);
      
      deriv = drdbt(f(j),bt(j,1:length(t)))*1E3;
      
      dr = (x/365).*all_b(j,2);
      lbt(:,j) = dr./deriv';
%      lbt(:,j) = (dr./(deriv(j)*1E3))./(x(end)/365);
   else
      lbt(:,j) = NaN;      
   end
end

full_bt_resid = all_bt_resid + lbt;
   

start_time = datenum(2002,9,1,0,0,0);
stop_time = datenum(2018,8,31,0,0,0);
dmtime = datenum(1958,1,1,0,0,rtime_mean);
dmtime = nanmean(dmtime,2);
ndi = find( dmtime >= start_time & dmtime <= stop_time);
nd = length(ndi);

ig = 1:2645;
for i=ig
   k = remove_6sigma(robs(:,i));
   j = remove_6sigma(robs(k,i));
n   ind(i).k = k(j);
end

for i=ig
   it = ind(i).k;
   fittime = dmtime(it);
   fity = squeeze(robs(it,i));
   x = fittime - fittime(1);
   y = squeeze(fity);
   dr = (x/365).*all_b(i,2);
   rnew(it,i) = robs(it,i) -dr;
end

