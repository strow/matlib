function dateout = mattime2tai(varargin)
% function dateout = mattime2tai(mattime)
% function dateout = mattime2tai(mattime, year0)
% function dateout = mattime2tai(yyyy,mm,dd,HH,MM,SS)
% function dateout = mattime2tai(yyyy,mm,dd,HH,MM,SS, year0)
%
% Convert MATLAB time to TAI time
%
% year0 - is the base year. AIRS - 1993, IASI/CrIS - 2000
%         If not present, will use 1993
%
% See also:  iasi2mattime, tai2mattime, mattime2iasi

% Written by Paul Schou - 10 Sep 2009  (paulschou.com)
%            Breno Imbiriba - 2013.05.10

if nargin == 1
  dateout = (varargin{1}-datenum(1993,1,1,0,0,0))*86400;
elseif nargin == 2
  dateout = (varargin{1}-datenum(varargin{2},1,1,0,0,0))*86400;
elseif nargin == 6
  dateout = (datenum(varargin{:})-datenum(1993,1,1,0,0,0))*86400;
elseif nargin == 7
  dateout = (datenum(varargin{1:6})-datenum(varargin{7},1,1,0,0,0))*86400;
end
