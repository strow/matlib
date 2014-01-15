function p=Prof_join_arr(parr)
% function p=Prof_join_arr(parr)
%
% Take an array of profiles parr(:) and transform into one profile structure p of arrays.
%
% See that parr(:) must also be a valid profile structure, with the last dimension equal in all fields.

% create final structure, with zeros.

nl=zeros(size(parr));
fnames1=fieldnames(parr(1));
nfields=length(fnames1);

% Do we have an empty structure?
if(nfields==0)
  p=struct();
  return
end
  
for i=1:length(nl);
  n1l(i)=length(getfield(parr(i),fnames1{1}));
end

%ntot=sum(nl);

p=struct();
for i=1:nfields
  %sz=size(getfield(parr(1),fnames1{i}));
  p=setfield(p,fnames1{i},[]); %zeros([sz(1:length(sz)-1) ntot])); 
end
% Copy the data:

%for iarr=1:length(parr) % loop over structures
%  for iff=1:length(fnames1) % loop over fields
%
%    % get size of field;
%    sz=size(getfield(parr(1),fnames1{iff})); 
%
%    % concatenate at the last dimension.
%    tfield=getfield(parr(iarr),fnames1{iff});
%    p=setfield(p,fnames1{iff},cat(length(sz),getfield(p,fnames1{iff}),tfield));
%  end
%end

for iff=1:length(fnames1) % loop over fields

  % transform each field into an array:

  % get the rank of the array.
  frank=numel(size(parr(1).(fnames1{iff})));

  % concat data at the last dimension (the rank)
  f0=cat(frank,parr(:).(fnames1{iff}));


  % add to the final structure;

  p.(fnames1{iff})=f0;
 
end


end

