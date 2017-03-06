function ht = p2hFAST(pin,airslevels,airsheights)
%this function takes pressure in mb, and converts to height in m

% load airsheights.dat
% load airslevels.dat

h = airsheights;
p = airslevels;
for ii = 1:100
  pavg(ii) = (p(ii+1)-p(ii))/log(p(ii+1)/p(ii));
end
if length(h) ~= length(pavg)
  fprintf(1,'length(height) length9(pavg) = %3i %3i \n',length(h),length(pavg))
  error('in p2hFAST length(h) ~= length(pavg)')
end  
%whos h pavg
%plot(h,pavg)

ht = interp1(pavg,h,pin);

if ((isnan(ht)) | ht > 7.05e4)
  ht  =  8.09e4;
end
