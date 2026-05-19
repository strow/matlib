function raaY = interp1_loop(xi,yi,raaX,method)

%% xi,yi are vectors
%% raaX is array, so raaY is array

[m1,n1] = size(xi);
[m2,n2] = size(yi);
if sum(abs(size(xi)-size(yi))) ~= 0
  error('interp1_loop.m : xi and yi are differenr sizes');
end

[m,n] = size(raaX);
if nargin == 3
  for ii = 1 : n
    raaY(:,ii) = interp1(xi,yi,raaX(:,ii));
  end
else
  for ii = 1 : n
    raaY(:,ii) = interp1(xi,yi,raaX(:,ii),method);
  end
end
