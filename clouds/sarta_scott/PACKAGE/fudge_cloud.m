function [ciout,cwout,fraci,fracw,fraciw] = fudge_cloud(cfrac,cc,ciwc,clwc);

% from ~hannon/fudge_cloud.m
% function [ciout,cwout,fraci,fracw,fraciw] = fudge_cloud(cfrac,cc,ciwc,clwc);
%
% Convert ECMWF cloud data into two fudge cloud profiles and cloud
% fractions.
%
% Input:
%    cfrac = [   1 x nobs] total cloud cover fraction
%    cc    = [nlev x nobs] level cloud cover fraction
%    ciwc  = [nlev x nobs] ice cloud amount profile
%    clwc  = [nlev x nobs] water cloud amount profile
%
% Output:
%    ciout = [nlev x nobs] ice cloud profile
%    cwout = [nlev x nobs] water cloud profile
%    fraci = [1    x nobs] ice cloud fraction
%    fracw = [1    x nobs] water cloud fraction
%    fraciw= [1    x nobs] ice+water cloud fraction
%

% Created: 28 Feb 2008, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nlev,nobs] = size(cc);


% Do ice cloud
casum = sum(ciwc); % total ice cloud amount
caf = ciwc .* cc;
fraci = sum( caf ) ./ casum;
ciout = caf ./ (ones(nlev,1)*fraci);


% Do water cloud
casum = sum(clwc); % total ice cloud amount
caf = clwc .* cc;
fracw = sum( caf ) ./ casum;
cwout = caf ./ (ones(nlev,1)*fracw);


% Fraction with both clouds (unknown so make a guess)
fraciw = min([fraci, fracw]);


%%% end of function %%%
