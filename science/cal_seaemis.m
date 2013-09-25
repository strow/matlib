function [nemis, efreq, seaemis] = cal_seaemis(zang);
%
% function [nemis, efreq, emis] = cal_seaemis(zang);
%
% Calculates the sea surface emissivity (for no wind) in the IR.
%
% Input:
%    zang = (1 x nang) desired zenith angles (0-70 degrees)
%
% Output:
%    nemis = (1 x nang) number of emissivity points
%    efreq = (nemis x nang) frequencies for emissivity points (cm^-1)
%    emis  = (nemis x nang) emissivity points
%

% Created: 9 October 2001, Scott Hannon
% Last updated: 1 February 2002, Scott Hannon: renamed zang (was scanang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check zang
sang=abs(zang);
ii=find(sang > 70);
sang(ii)=70;
nang=length(sang);
[irow,icol]=size(sang);
if (nang > 1 & irow == 1)
   sang=sang';  %  want sang to be (nang x 1) for interpolation
end

% Load sea surface emissivity data
load sea_emis_nowind
emis=emis';  % want emis to be (nang x nemis) for interpolation

% Convert angles to secants
secang=1./cos(angles*pi/180)';
secsang=1./cos(sang*pi/180);

% Interpolate emissivity data to desired angles
nemis=length(wn)*ones(1,nang);
efreq=wn*ones(1,nang);
seaemis=interp1( secang, emis, secsang, 'linear')'; % want (nemis x ang)

%%% end of function %%%
