function calcifp(idef)

% function calcifp(idef)
% 
% calcifp() computes values used by Fourier convolvers such as
% kifft.m, kifft.c, kfft.m, fconvkc.m, fconv.m, fconvkc.c, and
% fconv.c, as well as various test and display programs.
%
% INPUT (file)
%
% The input argument to calcifp is the name of a a short m-file
% containing instrument parameters; the argument idef should be
% the name of an m-file without the .m extension, e.g., "hisB1" 
% for "hisB1.m".
% 
% This m-file should set the following matlab variables:
% 
%   v1  - band start wavenumber
%   v2  - band end wavenumber
%   L1  - real optical path length
%   dvc - specified channel spacing, or 1/(2*L1)
%
% dvc can be set to something less than 1/2L, if desired, in which
% case a zero-filled path length L2 will be used to to interpolate
% the convolved channels.  The m-file <idef>.m can also reset various
% other values; see the default settings in the code, below.
% 
% OUTPUTS (globals)
%    
% Because of the large number of variables set by calcifp(),
% values are returned as globals; the calling program needs 
% only to declare global those values it will actually use.
% (Most routines need only a few of these values; however it
% is convenient to do all these calculations in one place.)
%
%   v1       band low wavenumber
%   v2       band high wavenumber
%   vmax     transform max wavenumber
%   rolloff  rolloff at band edges (1/cm)
% 
%   L1       true optical path length
%   L2       L2 > L1 for dvc < 1/(2*L1)
%   Lcut     path storage truncation point
%   Lmax     max optical path, 1/(2*dv)
% 
%   dvk      kcarta wavenumber spacing
%   dv       vlaser-scaled wavenumber spacing
%   dvc      convolved channel spacing
%   ddk      kcarta-scaled distance increment
%   dd       vlaser-scaled distance increment
% 
%   v1ind    v1 index (in vlaser point scaling)
%   v2ind    v2 index (in vlaser point scaling)
%   vmaxind  vmax index (in vlaser point scaling)
%   vLpts    v2ind-v1ind+1 (number of vlaser pts)
% 
%   L1ind    L1 index
%   L2ind    L2 index
%   Lcutind  Lcut index (last point saved)
%   Lmaxind  Lmax index (same as vmaxind)
% 
%   v1cind   v1 index (in channel point scaling)
%   v2cind   v2 index (in channel point scaling)
%   nchans   v2cind-v1cind+1 (number of channels)
% 
%   v1kind   v1 index (in kcarta point scaling)
%   v2kind   v2 index (in kcarta point scaling)
%   kcpts    v2kind-v1kind+1 (number of kcarta pts)

% H. Motteler, 8/24/98


global v1 v2 vmax rolloff
global L1 L2 Lcut Lmax
global dvk dv dvc ddk dd
global v1ind v2ind vmaxind vLpts
global L1ind L2ind Lcutind Lmaxind
global v1cind v2cind nchans
global v1kind v2kind kcpts

% defaults that can be changed in parameter m-file
%
npts   = 2^19;	% for an npts+1 cosine transform (via a 2*npts fft)
vlaser = 15799;	% HeNe laser frequency
vsf    = 4;	% vlaser scaling factor
Lcut   = 2.5;   % path length storage cutoff 
rolloff = 4;	% rolloff at band edge (1/cm)

% read user-defined parameter assignments
%
eval(idef)

% kcarta parameters
%
dvk = .0025;		% kcarta point spacing
k = nextpow2(v2/dvk);	% points for uninterpolated cosine transform
ddk = 1/(2^(k+1)*dvk);  % kcarta effective distance increment 
v1kind = round(v1/dvk) + 1;	% kcarta band low dvk index
v2kind = round(v2/dvk) + 1;	% kcarta band high dvk index
kcpts = v2kind - v1kind + 1;	% kcarta spaced points from v1 to v2

% vlaser and npts derived parameters
%
dd = vsf / vlaser;	% vlaser-scaled step size
dv = 1 / (2*npts*dd);   % vlaser-scaled wavenumber spacing

Lmax = 1 / (2*dv);	% max inverse transform path length
vmax = 1 / (2*dd);	% max inverse transform wavenumber

v1ind = round(v1/dv) + 1;	% band low index
v2ind = round(v2/dv) + 1;	% band high index
vmaxind = round(vmax/dv) + 1;	% band max index (should be npts+1)
vLpts = v2ind - v1ind + 1;	% vlaser spaced points from v1 to v2

% path length and channel parameters
%
L2 = 1/(2*dvc);			% L1 < L2 for "oversampled" dvc

L1ind = round(L1/dd) + 1;	% points in truncated interferogram
L2ind = round(L2/dd) + 1;	% points to L2
Lcutind = round(Lcut/dd) + 1;	% points to Lcut (last point saved)
Lmaxind = round(Lmax/dd) + 1;	% points to Lmax (same as vmaxind)

v1cind = round(v1/dvc) + 1;	% v1 channel index
v2cind = round(v2/dvc) + 1;	% v2 channel index
nchans = v2cind - v1cind + 1;	% number of channels

% adjust nominal v and L values
% 
v1 = (v1ind-1)*dv;
v2 = (v2ind-1)*dv;
L1 = (L1ind-1)*dd;
L2 = (L2ind-1)*dd;

% assorted sanity checks
%
if v2 > vmax
  fprintf(2, 'calcifp WARNING: v2 > vmax\n')
end

if dv > 4*dvk
  fprintf(2, 'calcifp WARNING: dv >> dvk\n')
end

if 2*dv <  dvk
  fprintf(2, 'calcifp WARNING: dv << dvk\n')
end

if vmaxind ~= npts + 1
  fprintf(2, 'calcifp ERROR: vmaxind=%d npts=%d\n', vmaxind, ntps);
  error('vmaxind should equal npts + 1');
end

