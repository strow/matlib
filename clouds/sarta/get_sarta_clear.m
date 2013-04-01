%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

klayers = run_sarta.klayers_code;
sarta   = run_sarta.sartaclear_code;

if ~exist(klayers,'file')
  error('klayers exec done not exist')
end
if ~exist(sarta,'file')
  error('sarta cloud exec done not exist')
end

fip = mktemp('temp.ip.rtp');
fop = mktemp('temp.op.rtp');
frp = mktemp('temp.rp.rtp');
ugh = mktemp('ugh');

rtpwrite(fip,h,ha,prof,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
  eval(klayerser);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh];
  eval(sartaer);
[headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh]; eval(rmer);
