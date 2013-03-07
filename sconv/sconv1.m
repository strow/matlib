
function [rout, fout] = sconv1(rin, fin, clist, mfile);

% NAME
%
%   sconv1 -- apply a sparse SRF convolution matrix
%
% SYNOPSIS
%
%   function [rout, fout] = sconv1(rin, fin, clist, mfile);
%
% INPUTS
%
%   rin    - m x n input data
%   fin    - m - vector of input frequencies
%   clist  - output channel list 
%   mfile  - sparse convolution .mat file
% 
% OUTPUTS
%
%   rout   - k x n convolved ouput data
%   fout   - k-vector of output frequencies
% 
% DESCRIPTION
% 
%   sconv1 calculates rout = Cmat(clist, flist) * rin, where
%   flist is the indices of fin in Cfin 
%
% H. Motteler, 15 Jan 02
%

% guarantee input vectors are columns
fin = fin(:);
clist = clist(:);

% check rin dimensions
[m,n] = size(rin);
if m ~= length(fin)
  error('number of rows of rin do not match length of fin')
end

% set a default matlab convolution file
if nargin < 4
  mfile = '/asl/data/airs/srf/srftablesV10.mat';
end

% load the convolution file, define the variables
%
%   Cmat    - sparse convolution matrix
%   Cfin    - Cmat input frequency list
%   Cfout   - Cmat output frequency list
%   Ciout   - Cmat output channel number
% 
eval(sprintf('load %s', mfile));

% find indices of fin in Cmat
rind = interp1(Cfin, 1:length(Cfin), fin, 'nearest');

% find indices of clist in Cmat
cind = interp1(Ciout, 1:length(Ciout), clist, 'nearest');

% do the convolution
rout = Cmat(cind, rind) * rin;

% return output frequency list
fout = Cfout(cind);

