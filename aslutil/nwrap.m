% FUNCTION NWRAP(X,N)
%
%  Adds columns on the edges from the reverse side, reversed by nchomp
%
%x =
%     1     2     3     4
%     5     6     7     8
%> nwrap(x)
%ans =
%     4     1     2     3     4     1
%     8     5     6     7     8     5
%
% N - determines the number of columns to wrap

% Written by Paul Schou (paulschou.com) - 3 June 2009
function y = nwrap(x,n)

if(nargin == 1)
    n = 1;
end
y = [x(:,end-(n-1):end) x x(:,1:n)];