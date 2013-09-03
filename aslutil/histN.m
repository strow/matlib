function histmat = histN(varargin)
% histN is a function to bin data over any given sets of dimensions.  This
%   function can be used to evaluate data that has variations in multiple
%   dimensions and make plots.
%
% histmat = histN(dataA, binsA)
% histmat = histN(dataA, dataB, binsA, binsB)
% histmat = histN(dataA, dataB, ... binsA, binsB ...)
% 
% Example:
%   [a]=histN(rand(10)*10,rand(10)*4,0:10,1:4)
%
% Note: A warning is given if data points have been left out of the selection.
%
% See also:  GSTATS

% Written by Paul Schou -- 10 April 2013

% checks for even number of inputs
if mod(nargin,2) ~= 0 || nargin == 0
  error('Number of inputs must be even:  histN(valA, valB, ... binA, binB ...)')
end

outsize = [];
for i = (nargin/2+1):nargin
  % checks for dataset uniformity
  if ~isequal(size(varargin{1}),size(varargin{i-nargin/2}))
    error(['Input dataset ' num2str(i-nargin/2) ' is not the same size as the first dataset'])
  end
  % checks for bin vectorizations
  if ~isvector(varargin{i}) || numel(varargin{i}) == 1
    error(['Bins for dimension ' num2str(i-nargin/2) ' is not a vector'])
  end
  outsize = [outsize length(varargin{i})-1];
end

% predeclare the bin set matrix
%histmat = zeros(outsize);
bins = zeros(numel(varargin{1}),length(outsize));
bin_zaps = zeros(1,length(outsize));

% the remaining part of this code is taken out of gstats and simplified
for i = 1:nargin/2
  [x bins(:,i)] = histc(varargin{i}(:),varargin{i+nargin/2});
  % remove any infinities or reasons for the N+1 bin
  bins(bins(:,i) >= length(varargin{i+nargin/2}),i) = 0;
  % keep track of how many points did not make a bin
  bin_zaps(i) = sum(bins(:,i) == 0);
end

% remove zeros
if any(bins(:) == 0)
  disp(['Warning: ' num2str(sum(bins(:)==0)) ' of ' num2str(size(bins,1)) ' datapoints not counted in the selected bin sizes.  Bin drops: {' num2str(bin_zaps,'%g ') '}'])
  bins = bins(~any(bins == 0,2),:);
end

% create the output matrix
histmat = accumarray(bins,1,outsize);
