function a = sinxx(x)

% function a = sinxx(x)
% 
% returns sin(x)/x, 1 at x = 0

d = x==0;
a = sin(x) ./ (x+d);
a = a .* (1-d) + d;

