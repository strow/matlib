function ht = p2hFAST(pin,airslevels,airslayers,airsheights)
%this function takes pressure in mb, and converts to height in m

global iWhichInterp  %% 0 = matlab interp1, 1 = interp1qr, set in driver_sarta_cloud_rtp.m

% load airsheights.dat
% load airslevels.dat

h = airsheights;
p = airslevels;
pavg = airslayers;

%%pavg = zeros(100,1);
%%for ii = 1:100
%%  pavg(ii) = (p(ii+1)-p(ii))/log(p(ii+1)/p(ii));
%%end
%pavgN = p(1:end-1)-p(2:end);
%pavgD = log(p(1:end-1)./p(2:end));
%pavg  = pavgN./pavgD;

%if length(h) ~= length(pavg)
%  fprintf(1,'length(height) length9(pavg) = %3i %3i \n',length(h),length(pavg))
%  error('in p2hFAST length(h) ~= length(pavg)')
%end  
%plot(h,pavg)

if iWhichInterp == 0
  ht = interp1(pavg,h,pin);
elseif iWhichInterp == 1
  [mm,nn] = size(pin);
  if mm*nn == 1
    ht = interp1qr(pavg,h,pin);
    %subplot(121); loglog(pin,ht0,'b',pin,ht,'r'); 
    %subplot(122); semilogy(ht0-ht,pin,'r'); 
    %pause(0.1)
  elseif mm*nn > 1
    ht = interp1qr(pavg,h*ones(1,nn),pin);
  end
end

%if ((isnan(ht)) | ht > 7.05e4)
%  ht  =  8.09e4;
%end

bad = find(isinf(ht) | isnan(ht) | ht > 8.09e4);
ht(bad) = 8.09e4;

