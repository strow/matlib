function [hjunk,hajunk,pjunk,pajunk] = get_sarta_cloud100layer_klayersONLY(h,ha,prof,pa,run_sarta);

%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

%% script is same as get_sarta_cloud except we call klayers100 and sarta100
%% calling klayers100 means we need to specify that we want gases 201,202 to be output

klayers = run_sarta.klayers_code;

if ~exist(klayers,'file')
  error('klayers exec done not exist')
end

fip = mktemp('temp.ip.rtp');
fop = mktemp('temp.op.rtp');
ugh1 = mktemp('ugh1');

gas_str = 'nwant=10 listg=1,2,3,4,5,6,9,12,201,202 ';

rtpwrite(fip,h,ha,prof,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' '  gas_str ' >& ' ugh1];
  eval(klayerser);
  [hjunk,hajunk,pjunk,pajunk] = rtpread(fop);

rmer = ['!/bin/rm ' fip ' ' fop ' ' ugh1 ]; eval(rmer);

for ix = 1 : length(prof.stemp)
  xnlevs = prof.nlevs(ix);
  xplevs = prof.plevs(1:xnlevs,ix);
  xcc    = prof.cc(1:xnlevs,ix);

  ynlevs = pjunk.nlevs(ix);
  yplevs = pjunk.plevs(1:ynlevs,ix);
  yplevsA = yplevs(1:end-1)-yplevs(2:end);
  yplevsB = log(yplevs(1:end-1)./yplevs(2:end));
  yplevs  = yplevsA./yplevsB;
  if iWhichInterp == 0
    ycc     = interp1(log(xplevs),xcc,log(yplevs),[],'extrap');
  elseif iWhichInterp == 1
    ycc     = interp1qr(log(xplevs),xcc,log(yplevs));
  end
  ycc(ycc > 1) = 1;
  ycc(ycc < 0) = 0;
  pjunk.cc(1:101,ix) = 0;
  pjunk.cc(1:ynlevs-1,ix) = ycc';
end

yplevsA = pjunk.plevs(1:100,:) - pjunk.plevs(2:101,:);
yplevsB = log(pjunk.plevs(1:100,:) ./ pjunk.plevs(2:101,:));
yplevs  = yplevsA./yplevsB;
ycc      = pjunk.cc(1:100,:);
ix = 1 : min(5,length(prof.stemp));
  plot(ycc(:,ix),yplevs(:,ix),'b',prof.cc(:,ix),prof.plevs(:,ix),'r')
  
clear xnlevs xplevs xcc ynlevs yplevs* ycc
