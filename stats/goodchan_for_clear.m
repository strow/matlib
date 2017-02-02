function ig = goodchan_for_clear(count);

[d,f] = size(count);
if f ~= 2378 
   ig = NaN;
   disp('Error:  count must be ndays x 2378')
   return
end

mc = nanmean(count);

ig = find( mc > 0.9*max(mc) );

% Get rid of wired channels, above tests won't catch them
ibad = [121 122 133 134];
ig = setdiff(ig,ibad);
