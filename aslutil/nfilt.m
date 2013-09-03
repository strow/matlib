% FUNCTION NFILT(DATA, MODULUS) - Version 1.2 (2009)
%
% This noise filter will look for single point noise events and remove them by 
% averaging the eight points closest to the noisy point like so:
%
% a=[ 1  1  1  1
%     1  15 1  1
%     1  1  1  1
%     1  1  1  1 ]
% will result in:
% nfilt(a) 
% ans=[ 1  1  1  1
%       1  1  1  1
%       1  1  1  1
%       1  1  1  1 ]
%
% MODULUS - used to work with data that loops over a peroid (monthly => 12)

% Written by Paul Schou (paulschou.com) - September 2006
%   March 2007 - updated by Paul Schou to speed operations up with vectorizing
%   3 June 2009 - updated by Paul Schou to handle data that loops over a 
%       period for example yearly data (modval = 12) 
function [out] = nfilt(indat,modval,ns)

if nargin < 2
  ns = 2;
end

dat = [indat(1,1) indat(1,:) indat(1,end);indat(:,1) indat indat(:,end);indat(end,1) indat(end,:) indat(end,end)];

if(nargin > 1 & ~isempty(modval))
    dat = mod(dat,modval);
end

% Find the mean of a pixel and the 8 points closest to it
m = (dat(1:end-2,1:end-2) + dat(1:end-2,2:end-1) + dat(1:end-2,3:end) + dat(2:end-1,3:end) + ...
	dat(3:end,3:end) + dat(3:end,2:end-1) + dat(3:end,1:end-2) + dat(2:end-1,1:end-2))/8;

% Find the standard deviation of these same points.
s = sqrt(((dat(1:end-2,1:end-2)-m).^2 + (dat(1:end-2,2:end-1)-m).^2 + (dat(1:end-2,3:end)-m).^2 + (dat(2:end-1,3:end)-m).^2 + ...
	(dat(3:end,3:end)-m).^2 + (dat(3:end,2:end-1)-m).^2 + (dat(3:end,1:end-2)-m).^2 + (dat(2:end-1,1:end-2)-m).^2 )/8);

if(nargin > 1)
    % save the results
    m1=m;s1=s;

    % now center the data around 180 degrees out of phase
    if(nargin > 1 & ~isempty(modval))
      dat = mod(dat - modval/2, modval) + modval/2;
    end

    % Find the mean of a pixel and the 8 points closest to it
    m = (dat(1:end-2,1:end-2) + dat(1:end-2,2:end-1) + dat(1:end-2,3:end) + dat(2:end-1,3:end) + ...
        dat(3:end,3:end) + dat(3:end,2:end-1) + dat(3:end,1:end-2) + dat(2:end-1,1:end-2))/8;

    % Find the standard deviation of these same points.
    s = sqrt(((dat(1:end-2,1:end-2)-m).^2 + (dat(1:end-2,2:end-1)-m).^2 + (dat(1:end-2,3:end)-m).^2 + (dat(2:end-1,3:end)-m).^2 + ...
        (dat(3:end,3:end)-m).^2 + (dat(3:end,2:end-1)-m).^2 + (dat(3:end,1:end-2)-m).^2 + (dat(2:end-1,1:end-2)-m).^2 )/8);

    sel = s1 < s;
    m(sel) = m1(sel);
    s(sel) = s1(sel);
end

% Subset the data for the part without edges
d = dat(2:end-1,2:end-1);

% select the points which are `noisy'
sel = find(abs(d-m) > ns * s | isnan(d));
% replace the selected points with the average of those around the point
d(sel) = m(sel);
% put the subsetted data back into the original data set
%dat(2:end-1,2:end-1) = d;

if(nargin > 1 & ~isempty(modval))
    % bring the span back to 0-modval
    d = mod(d, modval);
end

% return this new dataset
out = d;

