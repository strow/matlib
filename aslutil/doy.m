function ddd = doy(varargin)
% function ddd = doy(yyyy,mm,dd)
% function ddd = doy([yyyy mm dd ...])
% function ddd = doy(datenum)
%
% Return the day of the year
%
% Breno Imbiriba - 2013.05.06

if(nargin()==3)
  dn = datenum([varargin{1:3} 0 0 0]);
  d0 = datenum([varargin{1} 1 0 0 0 0]);
elseif(nargin()==1)
  if(numel(varargin{1})==1)
    dn = varargin{1};
    vv = datevec(dn);
    d0 = datenum([vv(1) 1 0 0 0 0]);
  elseif(numel(varargin{1})>=3)
    dn = datenum([varargin{1}(1:3) 0 0 0]);
    d0 = datenum([varargin{1}(1) 1 0 0 0 0]);
  else
    error('Wrong number of arguments');
  end
else
  error('Wrong number of arguments');
end

ddd = dn - d0;
end

