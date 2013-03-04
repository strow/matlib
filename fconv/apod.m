function a = apod(d, L1, atype, aparg);

% function a = apod(d, L1, atype, aparg);
%
% apply an apodization by name
%
% inputs
%   d     - distance; may be a vector
%   L1    - max path length
%   atype - apodization type (see below)
%   aparg - optional apodization parameter
%
% output
%   a - apodization of d
% 
% available apodizations
%   (atype name)	(description)
%   boxcar, box		boxcar (simple truncation at L)
%   triangle, tri	triange
%   hamming, ham	Hamming 
%   kaiser-bessel, kb	Kaiser-Bessel #aparg (aparg = 1 to 8)
%   norton-beer, nb	Norton-Beer #aparg (aparg = 1 to 8)
%   cosine, cos		cosine function (0 to pi/L, normalized to [0,1])
%   beer		Beer (from Wisconsin?)


if nargin == 2
  fprintf(2, 'apod() warning: using default, kaiser-bessel #6\n');
  atype = 'kaiser-bessel';
  aparg = 6;
end

switch atype

  case {'boxcar', 'box'},		a = boxapod(d, L1);
  case {'triangle', 'tri'},		a = triapod(d, L1);
  case {'hamming', 'ham'},		a = hamapod(d, L1);
  case {'kaiser-bessel', 'kb'},		a = kbapod(d, L1, aparg);
  case {'norton-beer', 'nb'},		a = nbapod(d, L1, aparg);
  case {'cosine', 'cos'},		a = cosapod(d, L1);
  case {'beer'},			a = beerapod(d, L1);
  case {'gauss'},			a = gaussapod(d, L1);
  otherwise, error(['unknown apodization: ', atype]);

end

