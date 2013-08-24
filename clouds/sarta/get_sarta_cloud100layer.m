%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

%% script is same as get_sarta_cloud except we call klayers100 and sarta100
%% calling klayers100 means we need to specify that we want gases 201,202 to be output

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

gas_str = 'nwant=10 listg=1,2,3,4,5,6,9,12,201,202 ';

rtpwrite(fip,h,ha,prof,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' '  gas_str ' >& ' ugh1];
  eval(klayerser);

sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh2];
  eval(sartaer);
try
  [headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
catch me
  me
  fprintf(1,'oops : error running sarta cloudy, look at error log %s \n',ugh2);
  %keyboard
  error('woof! try again!')
end
  
rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh1 ' ' ugh2]; eval(rmer);
