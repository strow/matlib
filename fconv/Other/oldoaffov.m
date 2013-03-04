function [oafreq,oasrf] = oaffov(fchan,fmin,fmax,df,opd,thetac,hfov,nslice);

%
%  function [oafreq,oasrf] = oaffov(fchan,fmin,fmax,df,opd,thetac,hfov,nslice);
%
%  Compute an interfermeter's spectral space SRF.  This routine handles
%  the case of an Off-Axis Finite Field Of View (OAFFOV).  The FFOV results
%  in both a broadening of the SRF, as well as a shift toward lower wave
%  numbers.  If thetac is not zero (ie the FFOV is off-axis), then the SRF
%  is not symmetric.  The basic sinc shape is due to the finite opd.
%  Note: this function needs to be used on a per channel basis.
%
%  Input:
%     fchan = channel center freq (cm-1)
%     fmin = min freq for output (cm-1)
%     fmax = max freq for output (cm-1)
%     df = freq spacing for output (cm-1)
%     opd = max OPD (cm)
%     thetac = the central off-axis angle (radians) >= 0
%     hfov = field of view half angle (radians) >= 0
%     nslice = number (integer>1) of slices for numeric integration
%
%     Note: both thetac and hfov should be small angles (maybe < 0.1).
%           nslice should be large (maybe > 1000).
%
%  Output:
%     oafreq = fmin to fmax with df pt spacing (cm-1)
%     oasrf = normalized (peak=1.0) response function
%

% The routine is based on info recieved from an e-mail from
% David Tobin on 3 June 1998.
% Created by Scott Hannon, November? 1998
% Last updated/corrected on 22 February 1999; fixed serious error in
%    integration step size and limits.

% The general formula for a zero FOV at some angle theta is given by:
%
%    sin(A)/A with A=2*pi*( freq - fchan*cos(theta) )*opd
%
% To get the full SRF, this has to be integrated over the circular
% FOV beam area.  Note the angle theta is the total off-axis angle
% of a point.
%
%    SRF = sum(i=1 to nslice) of { W_i * sin(A_i)/A_i }
%
% where W_i is half the length of the ith arc formed by the intersection
% of a circle of radius theta_i centered on the origin, and a circle of
% radius hfov centered on the x-axis at x=thetac
% Note: there is symmetry about the x-axis: the length of W_i above the
% x-axis is the same as the length below.  The total length is thus always
% twice the length above the x-axis, but this factor of two is a constant
% of the integration and drops out of the normalized response.
%
%    W_i = theta_i * phi_i, where phi is the azimuthal angle given by
%    phi_i = pi if abs( B_i ) > 1, otherwise phi_i =acos( B_i )
%    B_i = cos(phi_i) = (theta_i^2 - hfov^2 + thetac^2)/(2*theta_i*thetac)
%


% Assign freq points
oafreq=fmin:df:fmax;
npts=length(oafreq);

% Create work array
oasrf=zeros(1,npts);

% Assign theta step and initial theta value
mintheta=max(0, thetac - hfov);  % Min off-axis angle is zero
maxtheta=thetac + hfov;
dtheta=(maxtheta - mintheta)/(nslice - 1);
theta_i=mintheta;

% Loop over the slices
for i=1:nslice
   A_i=2*pi*( oafreq - fchan*cos(theta_i) )*opd;
   denom=2*theta_i*thetac;
   if (denom > 0)
      B_i = (theta_i^2 - hfov^2 + thetac^2)/denom;
      if ( abs(B_i) > 1 )
         phi_i=pi;
      else
         phi_i=acos(B_i);
      end
   else
      phi_i=pi;
   end
   W_i=abs(theta_i*phi_i); % The area of the slice = W_i >= zero
   slice=W_i*ones(1,npts);
   ind=find(A_i ~= 0);
   slice(ind)=W_i*sin(A_i(ind))./A_i(ind);
   oasrf=oasrf + slice;
   %
   % Update theta_i for next slice
   theta_i=theta_i + dtheta;

%   i,B_i,phi_i,W_i
%   plot(oafreq,slice)
%   pause

end

% Area normalize the output
oasrf=oasrf/sum(oasrf);
