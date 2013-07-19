% FUNCTION NCHOMP(X,N)
%
%  Does the reverse of nwrap and removes the added columns
%
%x =
%     4     1     2     3     4     1
%     8     5     6     7     8     5
%> nchomp(x)
%ans =
%     1     2     3     4
%     5     6     7     8
%
% N - determines the number of columns to wrap

% Written by Paul Schou (paulschou.com) - 3 June 2009
function y = nchomp(x,n)

if(nargin == 1)
    n = 1;
end
y = x(:,(n+1):end-n);