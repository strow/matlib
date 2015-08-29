%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

disp(' ')
fprintf(1,'  doing cloud sky SARTA calcs ....\n')
fprintf(1,'    sarta cloud = %s \n',run_sarta.sartacloud_code);

klayers = run_sarta.klayers_code;
sarta   = run_sarta.sartacloud_code;

if ~exist(klayers,'file')
  error('klayers exec done not exist')
end
if ~exist(sarta,'file')
  error('sarta cloud exec done not exist')
end

fip = mktemp('temp.ip.rtp');
fop = mktemp('temp.op.rtp');
frp = mktemp('temp.rp.rtp');
ugh1 = mktemp('ugh1');
ugh2 = mktemp('ugh2');

rtpwrite(fip,h,ha,prof,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh1];
  eval(klayerser);

sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh2];
  eval(sartaer);
try
  [headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
catch me
  me
  fprintf(1,'oops : error running sarta cloudy, look at ip/op rtp files %s %s \n',fip,fop)
  fprintf(1,'oops : error running sarta cloudy, look at error log %s \n',ugh2);
  %keyboard
  error('woof! try again!')
end
  
rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh1 ' ' ugh2]; eval(rmer);
