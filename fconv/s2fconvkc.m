function [rch, vch] = s2fconvkc(rkc, ifp, atype, aparg)

% function [rch, vch] = s2fconvkc(rkc, ifp, atype, aparg)
% 
% Fourier convolution of monochromatic/kCARTA spectra.  The
% convolution parameters are set by the specified "ifp" file.
% Version2 with npts =2^(pow2min+1)
%
% Inputs
%   rkc   - [kcpts x nprof] monochromatic spectra, where kcpts
%      matches matches the length of v1:dvk:v2 specified by ifp.
%   ifp   - [string] name of the interferometric parameter file
%     (without the .m, eg "tansotir"). At a minimum, variables
%     {L1, v1, v2} must be specified in the ifp file.
%   atype - [string] optional apodization type (default 'box').
%      Available choices: box, gauss, ham, beer, tri, cos, nb, kb
%   aparg - [1 x 1] optional apodization parameter for kb & nb
%
% Outputs
%   rch   - [nchans x nprof] convolved spectra
%   vch   - [nchans x 1] channel center wavenumbers
%

% Created: 25 Aug 1998, Howard Motteler
% Update: 08 Jan 2010, Scott Hannon - re-write to simplify and improve
%    robustness/accuracy for arbitraty dvc.  The original fconvkc can be
%    innaccurate if dvc is unrelated to dvk; this code is always OK.
%    The arguments are unchanged except for vch, which was previously
%    called "wch" and [1 x n].  This code is compatible with old ifp files,
%    but I recommend they be updated to specify vcmin & vcmax (these will
%    have no effect when used with the original fconvkc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check number of input arguments
if (nargin < 2)
   error('insufficient input arguments for xfconvkc')
end
if (nargin == 3)
   aparg = -1;
elseif (nargin == 2)
   atype = 'boxcar';
   aparg = -1;
end  


% Load interferometric parameters
clear global v1 v2
clear global L1
clear global vcmin vcmax
clear global dvk dvc
clear global rolloff
clear global kcpts
clear global nchans
%
scalcifp(ifp);
%
global v1 v2           % min & max kCARTA freq points {cm^-1}
global L1              % true optical path difference {cm}
global vcmin vcmax     % min & max output channel freqs {cm^-1}
global dvk dvc         % kCARTA and output freq point spacing {cm^-1}
global rolloff         % spectral rolloff distance {cm^-1}
global kcpts           % number of kCARTA points, ie length(v1:dvk:v2)
global nchans          % number of convolved channels


% Check supplied spectra
[m,nprof] = size(rkc);
if (m ~= kcpts)
   m
   kcpts
   error('number of spectral points in rkc does not match ifp file')
end


% Input kCARTA freq points
vk = v1:dvk:v2;


% Frequency points for FFT
dv = dvk;    % Stick with the kCARTA grid
%
v1ind = round(v1/dv + 1);
v2ind = round(v2/dv + 1);
%
v = ((v1ind:v2ind) - 1)*dv;


% Assign npts and nfft
pow2min = nextpow2(v2ind); % min required power of 2
pow2 = pow2min + 1;        % power of 2 to be used; pow2> = pow2min
npts = 2^pow2;             % number of points in zero padded mono spectrum
nfft = 2*npts;             % number of FFT points is always 2*npts


% Compute spectral rolloff filter for band edges
nro = round(rolloff/dv);  % number of non-unity rolloff points
nrom1 = nro - 1;
f0 = ( 1 + cos(pi + pi*(0:(nrom1))/nro) )/2;
rofilt = ones(nfft,1);
rofilt(v1ind:(v1ind+nrom1)) = f0;
rofilt((v2ind-nrom1):v2ind) = fliplr(f0);


% Declare output arrays
rch = zeros(nchans, nprof);
vch = (vcmin:dvc:vcmax)'; %'


% OPD points for left side of intf1
d = 1/(dv*2)*linspace(0,1,npts+1); % Note d is length nfft/2 + 1

% Freq grid for rall; only the first/left half is used/correct.
vnfft = (0:(nfft-1))*dv;


ind0l= find(d > L1);
% Note: last point for max(d) on the left is not repeated on the right
ii = length(ind0l) - 1;
ind0r = nfft + 2 - ind0l(1:ii);
ind0 = union(ind0l,ind0r);



% Build the apodization vector
apl = apod(d, L1, atype, aparg)'; %' left side apod vector, npts+1
apr = flipud(apl(2:npts));        % right side apod vector, npts-1
apvec = [apl; apr];                 % complete apod vector, 2*npts


% Loop over the profiles and do convolution
rv = zeros(nfft,1);  % work vector for current spectra 0:dv:vmax-dv
for ir =1:nprof
   if (dv == dvk)
      rv(v1ind:v2ind) = rkc(:,ir);
   else
      % Interpolate from kCARTA freq grid to FFT freq grid
      rv(v1ind:v2ind) = interp1(vk,rkc(:,ir),v,'linear');
   end

   % Apply spectral rolloff filter
   rv = rv .* rofilt;


   % Compute interferograsm
   % Note: rv has a correct spectrum on the left side, but the reverved
   % spectrum for the right side has been omitted (not strictly needed)
   intf = real( ifft(rv,nfft) );
   % "intf" is a 2 sided interferogram with ZPD at the ends and max OPD
   % in the middle.  Points at index 1 and nfft/2 + 1 are NOT repeated;
   % all other points on the left side are repeated on the right side
   % according to intf(2:(nfft/2)) = intf(nfft:-1:(nfft/2+2)).


   % Apply apodization
   intf = intf .* apvec;


   % Transform back into spectra (still at dv point spacing)
   rx = real( fft(intf,nfft) )*2;
   % Note: the 2 above is needed to account for the omitted reversed
   % spectrum on the right side of rv.


   % Interpolate to desired output channels
   % Note: dv spacing is very fine so linear interp should be adequate.
   rch(:,ir) = interp1(vnfft,rx,vch,'linear');

end

%%% end of function %%%
