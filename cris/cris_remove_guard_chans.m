function [fg,f,k] = cris_remove_guard_chans(restype);

% Input: restype = 'hires' or 'lowres'
% Ouptup: f (fg) = frequencies w/o (w) guard channels, k = indices to remove guard channels

addpath /asl/packages/ccast/motmsc/utils/
addpath /asl/packages/ccast/motmsc
addpath /asl/packages/ccast/motmsc/time


% Get proper frequencies for these data
switch restype
  case 'lowres'
    [n1,n2,n3,userLW,userMW,userSW, ichan] = cris_lowres_chans();
  case 'hires'
    [n1,n2,n3,userLW,userMW,userSW, ichan] = cris_highres_chans();     
  otherwise
    disp(['Bad restype;  Must either be ''hires'' or ''lowres'''])
end

fg = cris_vchan(2, userLW, userMW, userSW);
nchan = n1+n2+n3;

% Subset to real channels
k = find(ichan <= nchan);
f = fg(k);
f = f(:);
fg = fg(:);
k = k(:);

