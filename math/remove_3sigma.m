function [i] = remove_3sigma(x);

xstd  = nanstd(x(:));
xmean = nanmedian(x(:));

i = find( abs(x(:)-xmean) < 3*xstd);

%disp('       % Removed ');
%disp(100*(length(x(:)) - length(i))./length(x(:)));

return




