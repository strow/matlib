function s = apname(atype, aparg)

% function s = apname(atype, aparg)
%
% inputs
%  atype - apodization type
%  aparg - optional parameter
%
% output
%  s - fancy apodization name

switch atype

  case {'boxcar', 'box'},	s = 'boxcar';
  case {'triangle', 'tri'},	s = 'triangle';
  case {'hamming', 'ham'},	s = 'Hamming';
  case {'kaiser-bessel', 'kb'},	s = ['Kaiser-Bessel #', num2str(aparg)];
  case {'norton-beer', 'nb'},	s = ['Norton-Beer #', num2str(aparg)];
  case {'cosine', 'cos'},	s = 'cosine';
  case {'beer'},		s = 'Beer';
  case {'gauss'},		s = 'Gaussian';
  otherwise, error(['unknown apodization: ', atype]);

end

