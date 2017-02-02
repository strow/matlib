function ig = goodchan_onlybycount(count,cutoff);

[d,f] = size(count);
if f ~= 2378 
   ig = NaN;
   disp('Error:  count must be ndays x 2378')
   return
end

stdc = nanstd(count)-nanstd(count(:,1));
% Note most non-PC arrays stdc will be -0.04 

ig = find( stdc > -0.2 & stdc < cutoff );

% Get rid of wired channels, above tests won't catch them
ibad = [121 122 133 134];
ig = setdiff(ig,ibad);
