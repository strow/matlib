function [h12] = p2hFAST2(pin1,pin2,airslevels,airslayers,airsheights)
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
  h12 = interp1(pavg,h,[pin1 pin2]);
elseif iWhichInterp == 1
  h12 = interp1qr(pavg,h*ones(1,2),[pin1 pin2]');
end

%if ((isnan(h12)) | h12 > 7.05e4)
%  h12  =  8.09e4;
%end

bad = find(isinf(h12) | isnan(h12) | h12 > 8.09e4);
h12(bad) = 8.09e4;

