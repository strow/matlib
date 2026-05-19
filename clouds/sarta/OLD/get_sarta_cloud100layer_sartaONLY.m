%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

%% script is same as get_sarta_cloud except we call klayers100 and sarta100
%% calling klayers100 means we need to specify that we want gases 201,202 to be output

sarta   = run_sarta.sartacloud_code;

if ~exist(sarta,'file')
  error('sarta cloud exec done not exist')
end

fop = mktemp('temp.op.rtp');
frp = mktemp('temp.rp.rtp');
ugh2 = mktemp('ugh2');

rtpwrite(fop,hX,ha,profX,pa);
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
  
rmer = ['!/bin/rm ' fop ' ' frp ' ' ugh2]; eval(rmer);
