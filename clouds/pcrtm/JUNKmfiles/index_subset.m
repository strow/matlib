function p0 = index_subset(NN,pIN);

%  inds = 1 : 51;
%  inds = 1 : 101;

p0 = pIN;

inds = NN;

p0.efreq = p0.efreq(:,inds); 
p0.emis  = p0.emis(:,inds); 
p0.rho   = p0.rho(:,inds); 
p0.cpsize  = p0.cpsize(:,inds); 
p0.cpsize2 = p0.cpsize2(:,inds); 
if isfield(p0,'cfrac')
  p0.cfrac   = p0.cfrac(:,inds); 
end
if isfield(p0,'cfrac2')
  p0.cfrac2  = p0.cfrac2(:,inds); 
end
if isfield(p0,'cfrac12')
  p0.cfrac12 = p0.cfrac12(:,inds); 
end

p0.rtime   = p0.rtime(:,inds); 
p0.atrack  = p0.atrack(:,inds); 
p0.xtrack  = p0.xtrack(:,inds); 

if isfield(p0,'udef')
  p0.udef = p0.udef(:,inds); 
end

p0.rcalc = p0.rcalc(:,inds); 

