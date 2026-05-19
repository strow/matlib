function [prof,hlayers,players] = get_sarta_cloud(h,ha,prof0,pa,run_sarta);

tic
prof = prof0;

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

%fip = mktempS('temp.ip.rtp');
%fop = mktempS('temp.op.rtp');
%frp = mktempS('temp.rp.rtp');
%ugh = mktempS('ugh');
%ugh2 = mktempS('ugh2');

fip = mktemp('temp.ip.rtp');
fop = mktemp('temp.op.rtp');
frp = mktemp('temp.rp.rtp');
ugh = mktemp('ugh');
ugh2 = mktemp('ugh2');

%{
%% sometimes the files get toooooo large
%%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if isfield(prof,'rcalc')
  prof = rmfield(prof,'rcalc');
end
if isfield(prof,'rcalc_std')
  prof = rmfield(prof,'rcalc_std');
end
if isfield(prof,'rad_allsky_std')
  prof = rmfield(prof,'rad_allsky_std');
end
if isfield(prof,'rad_clrsky')
  prof = rmfield(prof,'rad_clrsky');
end
if isfield(prof,'sarta_cloud')
  prof = rmfield(prof,'sarta_cloud');
end
if isfield(prof,'sarta_rclear')
  prof = rmfield(prof,'sarta_rclear');
end
if isfield(prof,'sarta_xclear')
  prof = rmfield(prof,'sarta_xclear');
end
h
ha
pa
prof
fprintf(1,'removed pcrtm \n')
%%%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%}

rtpwrite(fip,h,ha,prof,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
  eval(klayerser);

iBad = -1;
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh2];
  eval(sartaer);
try
  [headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
catch me
  iBad = +1;
  me
  disp('in get_sarta_cloud.m')  
  which mktemp
  which get_sarta_cloud
  fprintf(1,'oops : error running sarta cloudy, look at ip/op rtp files %s %s \n',fip,fop)
  fprintf(1,'oops : error running sarta cloudy, look at error log %s \n',ugh2);
  lser   = ['!ls -lth ' fip ' ' fop]; eval(lser);
  catter = ['!more ' ugh2];  eval(catter)
  %keyboard
  error('woof! try again!')
end
  
if iBad < 0
  rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh ' ' ugh2]; eval(rmer);
end

hlayers = headRX2;
players = profRX2;
toc
prof.rcalc = profRX2.rcalc;

