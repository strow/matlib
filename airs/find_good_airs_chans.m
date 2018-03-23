% find_good_airs_chans.m
%
%   Filter AIRS channels for "good" channels, meaning good over 14
%   years.  Use random stats files to examine counts per channel.
%   When the rtp and stats are created the calflag is read and the
%   observation is ignored and no count is added.  Then, commented out
%   below I look by hand at the clear scene stats for any other
%   outliers.  This produces indices "kg" which are the supposedly
%   good channels.  Additional indices are created for good channels
%   with no ab state changes, and subsets of that set of data for the
%   value of a and b (0,1,2).

%  Main selection is "percent_count_needed".  That is coded up here as
%  the percent of max count per year.  The is commented out code to do
%  this per the whole mission.  If you do 90% for the whole mission,
%  you get back channels like 791.7, CO2 Q-branch, which was off for a
%  year+.  But, since it's ab state changed, it is not included in the
%  later indices (kg_fixed_ab).

addpath /asl/matlib/time
addpath /asl/matlib/airs

% Use random time series to get statistics on how often a channel is good
% Use tropical 
fstart = '/asl/s1/strow/Data/Work/Airs/Random/Data/Desc/statlat';
count = NaN(40,5058,2378);
for latid = 10:30
   latid
   g = load([fstart int2str(latid)],'count');
   count(latid,:,:) = g.count;
end
% Get a time
load([fstart int2str(20)],'rtime');
mtime = tai2dtime(rtime);

% Average over middle lats, slight sampling problems at the poles
count = squeeze(nanmean(count));

% Want 90% for each year, I did a version for 0.9 over all years (good_chans_2016_badyears)
percent_count_needed = 0.9;
myears = year(mtime);
for i = min(myears):max(myears)
   k = find( myears == i );
   mycount = nanmean(count(k,:));
   nig = find( mycount > percent_count_needed*max(mycount) );
   if i > min(myears)
      nig = intersect(ig,nig);
   end
   ig = nig;
end

% % An alternative, to give you more channels, like the 791 Q-branch,
% % this lets a whole year 1/14 pass with no data
% mycount = nanmean(count);
% nig = find( mycount > percent_count_needed*max(mycount) );
% ig = nig;

% Remove channels wired together
ibad = [121 122 133 134];
ig = setdiff(ig,ibad);

%% Now look in Stability stat files for bad channels
%% fstart = '/asl/s1/strow/Data/Work/Airs/Random/Data/Desc_fits/fit_bias_lat';
% fstart = '/asl/s1/strow/Data/Work/Airs/Random/Data/Desc_fits/fit_nucal_bias_lat';
% for latid = 1:40
%    latid
%    g = load([fstart int2str(latid)]);
%    dbt_bias(latid,:) = g.dbt;
%    dbt_bias_err(latid,:) = g.dbt_err;
%    lag_bias(latid,:) = g.lag;
%    bias_mean(latid,:) = g.all_b(:,1);
% end
% % plot(ig,nanmean(dbt_bias(:,ig)),'+-'); % Find outliers
% %
% % I found two: 497 and 1355
% % Looking at counts versus channel 1:
% %    497: randomly distributed fewer counts by 6 or so out of 250
% %   1355: off on Oct 9, 2011, from 2012:2014 6-10 fewer counts, OK in 2014-mid 2015, starting up 
% %           again in late 2015 with 6-10 fewer counts

% Results from commented out code above
ibad = [ 312 497 1355 ];  % Only seem to be in nucal stats
% Channel 312 went bad March 2014, still working kinda, but -1K offset then
ig = setdiff(ig,ibad);

% AB changes
%
% First when looking for std(ab) == 0 I find 21 channels that were
% changed for 1 day (09-Oct-2011).  I will consider them "stable"

% I bunch (11 days) of times are NaT.  Not sure why, always ignored,
% thus the if statement below
for i=1:length(mtime);
   x = get_ab_state(mtime(i));
   if length(x) == 2378
      ab(i,:) = x;
   end
end
ab = single(ab);

% I find 13 channels where ab changed for one day, Oct 9, 2011
% You can find these with 
%   k0 = find(nanstd(ab) > 0.02 & nanstd(ab) < 0.04);
% So, I define no ab changes with nanstd(ab) < 0.05;
k = find( nanstd(ab) < 0.05 );

% Rename ig to kg
kg = ig;  % Good channels, some may change ab state over time
kg_fixed_ab = intersect(k,ig); %Good channels with no ab state changes

% But, could still have some bad channels with ab > 2.  (One found below.)
k = find(nanmean(ab(:,kg_fixed_ab)) > 2);
kg_fixed_ab = setdiff(kg_fixed_ab,kg_fixed_ab(k));

% Now find ab = 0,1,2 within kg_fixed_ab
mab = nanmean(ab);
% The 5E-4 prevents the nanstd(ab) 0.02 to 0.04 being caught again
kab0 = find(mab(kg_fixed_ab) <= 5E-4); kg_ab0 = kg_fixed_ab(kab0);
kab1 = find(mab(kg_fixed_ab) == 1);    kg_ab1 = kg_fixed_ab(kab1);
kab2 = find(mab(kg_fixed_ab) == 2);    kg_ab2 = kg_fixed_ab(kab2);

% % Channel summary
% kg              good chans, ab can change
% kg_fixed_ab     good chans, no ab chanes
% kg_ab0          good chans, ab=0 (a+b)
% kg_ab1          good chans, ab=1 (a only)
% kg_ab2          good chans, ab=2 (b only)

% Note "newline" command is new in 2017a
note = ['kg = good chans' newline 'kg_fixed_ab = good chans with fixed ab' newline 'kg_ab0/1/2 = good chans with fixed ab of 0/1/2'];
note2 = ['Can relax good chans by changing percent_count_needed = 0.9,' newline 'or my doing a full mission percent count (code commented out)'];

% save /asl/matlib/airs/good_chans_2016  kg* note note2 

