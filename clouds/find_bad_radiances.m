function [good,bad] = find_bad_radiances(freqs,rads)

addpath /asl/matlib/aslutil

%% sometimes robs1 for some chans is -9999 or close to 0, or
%% sometimes pcrtm calcs or sarta calcs fail
%%
%% this routine find the fovs where the radiances are bad
%%
%% eg bad = find_bad_radiances(h.vchan,p.rad_clrsky);

tp = rad2bt(freqs,rads);
bad = find(abs(imag(tp)) > eps); oof = zeros(size(tp)); oof(bad) = 1; oof = sum(oof); bad = find(oof > 0);

[mm,nn] = size(rads);
good = 1 : nn;
good = setdiff(good,bad);

if length(bad) == 1
  plot(freqs,rads(:,bad))
elseif length(bad) > 1
  plot(freqs,rads(:,bad),'b'); hold on;   plot(freqs,nanmean(rads(:,bad)'),'r','linewidth',2); hold off
end

fprintf(1,'found %5i good and %5i bad rads \n',length(good),length(bad))

