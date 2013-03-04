function scalcifp(ifp)

% function scalcifp(ifp)
% 
% Computes values used by Fourier convolver fconvkc
%
% Input:
%    ifp - [string] name of a a short m-file containing instrument
%       parameters (eg "tansotir").  At a minimum, the m-file must set
%       the following variables:
%          v1  - kCARTA start wavenumber
%          v2  - kCARTA end wavenumber
%          L1  - true maximum optical path difference
%
% Output: (globals)
%    
% Because of the large number of variables set by scalcifp(),
% values are returned as globals; the calling program needs 
% only to declare global those values it will actually use.
%
%   v1       band low wavenumber {cm^-1}
%   v2       band high wavenumber {cm^-1}
%   vmax     transform max wavenumber {cm^-1}
%   rolloff  rolloff at band edges {cm^-1}
%   dvk      kCARTA wavenumber spacing {cm^-1}
%   dvc      convolved channel spacing {cm^-1}
%   vcmin    min convolved channel wavenumber {cm^-1}
%   vcmax    max convolved channel wavenumber {cm^-1}
% 
%   L1       true Optical Path Difference (OPD) {cm}
%   L2       OPD implied by dvc {cm}
%   Lk       OPD implied by dvk {cm}
%   dd       dvk implied distance increment {cm}
%
%   v1ind    v1 index (in 0:dvk:v1)
%   v2ind    v2 index (in 0:dvk:v2)
%   L1ind    L1 index (in 0:dd:L1)
%   L2ind    L2 index (in 0:dd:L2)
%
%   nchans   number of convolved channels
%   kcpts    number of kcarta points
%   NFFT     number of FFT points
%

% HISTORY:
% Created: 24 Aug 1998, Howard Motteler (calcifp.m)
% Update: 06 Jan 2010, Scott Hannon - complete re-write to simplify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare globals
% First clear them in case they already exist
clear v1 v2 vmax rolloff dvk dvc vcmin vcmax
clear L1 L2 Lk dd v2inv v2ind L1ind L2ind
clear nchans kcpts nfft
%
global v1 v2
global vmax
global rolloff
global dvk dvc
global vcmin vcmax
global L1 L2 Lk dd
global v1ind v2ind
global L1ind L2ind
global nchans kcpts npts

% Read user-defined parameter assignments
%
eval(ifp)
%
if (min(size(v1)) == 0)
   error('ifp file is missing required variable v1')
end
if (min(size(v2)) == 0)
   error('ifp file is missing required variable v2')
end
if (min(size(L1)) == 0)
   error('ifp file is missing required variable L1')
end


% Set defaults as needed
if (min(size(L2)) == 0)
   L2 = L1;
end
if (min(size(dvk)) == 0)
   dvk = 0.0025;  % kCARTA point spacing for the usual 605-2830 cm^-1 database
end
if (min(size(rolloff)) == 0)
   rolloff = 4;
end
if (min(size(dvc)) == 0)
   dvc = 1/(2*L2);
end
if (min(size(vcmin)) == 0)
   vcmin = dvc*ceil( (v1+rolloff)/dvc );
end
if (min(size(vcmax)) == 0)
   vcmax = dvc*floor( (v2-rolloff)/dvc );
end


% Sanity checks
if (v1 <= 0)
   error('v1 must be positive')
end
if (v2 <= 0)
   error('v2 must be positive')
end
if (v2 < v1)
   error('Must have v2 > v1')
end
if (dvk <= 0)
   error('dvk must be positive')
end
if (L1 <= 0)
   error('L1 must be positive')
end
if (L2 < L1)
   error('L2 must be >= L1')
end
if (rolloff < 0)
   error('rolloff may not be negative')
end
if (rolloff < 4)
   disp('WARNING: rolloff is smaller than recommended minimum=4')
end
if (vcmin < v1+rolloff)
   disp('WARNING: vcmin smaller than recommended minimum=v1+rolloff')
end
if (vcmax > v2-rolloff)
   disp('WARNING: vcmax larger than recommended maximum=v2-rolloff')
end


% Set computed variables
v1ind = round(v1/dvk + 1);
v2ind = round(v2/dvk + 1);
nchans = length(vcmin:dvc:vcmax);
kcpts = v2ind - v1ind + 1;
%
pow2 = nextpow2(v2ind);
if (pow2 > 22)
   disp('WARNING: very large number of points for FFT; calc will be slow!')
end
npts = 2^pow2;   % nfft = 2*npts
%
Lk = 1/(dvk*2);  % max inverse transform optical path difference
dd = Lk/npts;    % optical path difference point spacing
vmax = 1/(2*dd); % max inverse transform wavenumber
%
L1ind = round(L1/dd + 1);
L2ind = round(L2/dd + 1);


%%% end of function %%%
