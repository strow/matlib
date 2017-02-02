function ig = goodchan_stats(count,bias_std,maxstd);

%keyboard

count = count(:);
bias_std = bias_std(:);

k = find( count < 2);
count(k) = NaN;
i = remove_3sigma(count);
ib = setxor(i,1:2378);
count(ib) = NaN;

mcount = nanmean(count(280:1720));
% Keep channels within 0.5% of max
k = find( count/mcount > 0.9 );
% Get rid of 2 sets of channels wired together
ibad = [121 122 133 134];
ig = setdiff(k,ibad);

k = find( bias_std(ig) < maxstd );
ig = ig(k);
