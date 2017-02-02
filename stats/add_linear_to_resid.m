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
   