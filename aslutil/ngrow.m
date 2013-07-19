% FUNCTION dat = NGROW(dat)
%
% Eliminates holes in data where NaNs reside by replacing NaN pixels by the
% mean of the surrounding pixels.  3 or more surrounding NaN pixels are
% required for the replacement to take place.
%
% f =
%   NaN   NaN   NaN   NaN
%   NaN     1     4     7
%   NaN   NaN     1     2
%> ngrow(f)
%ans =
%   NaN   NaN     4   NaN
%   NaN     1     4     7
%   NaN     2     1     2

% Written by Paul Schou (paulschou.com) - 3 June 2009
function y = ngrow(x,count)
y=x;

if nargin < 2
  count = 2;
end

for i = 1:size(x,1)
    for j = 1:size(x,2)
        if(isnan(x(i,j)))
            t = x(max(i-1,1):min(i+1,end),max(j-1,1):min(j+1,end));
            if(sum(~isnan(t(:))) > count)
                y(i,j) = mean(t(~isnan(t(:))));
            end
        end
    end
end
