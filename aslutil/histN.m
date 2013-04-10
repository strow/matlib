function histmat = histN(varargin)

% checks for even number of inputs
if mod(nargin,2) ~= 0
  error('Number of inputs must be even, (valA, valB, ... binA, binB ...)')
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

for i = 1:nargin/2
  [x bins(:,i)] = histc(varargin{i}(:),varargin{i+nargin/2});
end

% remove zeros
if any(bins(:) == 0)
  disp(['Warning: ' num2str(sum(bins(:)==0)) ' of ' num2str(size(bins,1)) ' datapoints not counted in the selected bin sizes'])
  bins = bins(~any(bins == 0,2),:);
end

% create the output matrix
histmat = accumarray(bins,1,outsize);
