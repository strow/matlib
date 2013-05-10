function dateout = tai2mattime(datein, y0)
%function dateout = tai2mattime(datein, y0)
%  TAI to mattime conversion
%
%  datein - TAI time from year0 (seconds past y0/1/1 0:0:0)
%  year0  - start year (if not present, will use 1993)
%           Use AIRS - 1993, IASI/CrIS - 2000
%  
% See also:  iasi2mattime, mattime2tai, mattime2iasi

% Paul Schou/Breno Imbiriba - 2013.05.10

  if(nargin==1)
    y0=1993;
  end

  dateout = datenum(y0,1,1,0,0,double(datein));

end
