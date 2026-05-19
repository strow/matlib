function [profOUT,hlayers,players] = get_sarta_clear(h,ha,prof0,pa,run_sarta)

tic

%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

%%%%%%%%%%%%%%%%%%%%%%%%%
profOUT = prof0;

%%%%%%%%%%%%%%%%%%%%%%%%%
prof = prof0;
  disp('explicitly setting cloud fields in working copy of prof to 0 and -9999')
  prof.cngwat = ones(size(prof.stemp)) * 0;
  prof.cfrac  = ones(size(prof.stemp)) * 0;
  prof.cpsize = ones(size(prof.stemp)) * -9999;
  prof.cprtop = ones(size(prof.stemp)) * -9999;
  prof.cprbot = ones(size(prof.stemp)) * -9999;
  prof.ctype  = ones(size(prof.stemp)) * -9999;

  prof.cngwat2 = ones(size(prof.stemp)) * 0;
  prof.cfrac2  = ones(size(prof.stemp)) * 0;
  prof.cpsize2 = ones(size(prof.stemp)) * -9999;
  prof.cprtop2 = ones(size(prof.stemp)) * -9999;
  prof.cprbot2 = ones(size(prof.stemp)) * -9999;
  prof.ctype2  = ones(size(prof.stemp)) * -9999;

  prof.cfrac12  = ones(size(prof.stemp)) * 0;
%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
fprintf(1,'  doing clear sky SARTA calcs ....\n')
fprintf(1,'    sarta clear = %s \n',run_sarta.sartaclear_code);

klayers = run_sarta.klayers_code;
sarta   = run_sarta.sartaclear_code;

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

%disp('inside /umbc/rs/pi_sergio/WorkDirDec2025/sergio_matlib/matlib/clouds/sarta/get_sarta_clear.m calling mktemp')
fip = mktemp('temp.ip.rtp');
fop = mktemp('temp.op.rtp');
frp = mktemp('temp.rp.rtp');
ugh = mktemp('ugh');
ugh2 = mktemp('ugh2');

%disp(' ')
%disp('inside /umbc/rs/pi_sergio/WorkDirDec2025/sergio_matlib/matlib/clouds/sarta/get_sarta_clear.m')
%which mktemp
%fprintf(1,'\n fip = %s \n ugh = %s \n',fip,ugh);

rtpwrite(fip,h,ha,prof,pa);
%printarray([min(prof.rlon) max(prof.rlon) min(prof.rlat) max(prof.rlat)],'in get_sarta_clear.m : min/max rlon  min.max rlat')
%'nanabooboo in /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/matlib/clouds/sarta/get_sarta_clear.m'
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
  eval(klayerser);
try
  [headRX2 hattrR2 profRX2 pattrR2] = rtpread(fop);
catch me
  me
  fprintf(1,'oops : error running klayers, look at error log %s \n',ugh);
  %keyboard
  error('woof! try again!')
end

sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh2];
  eval(sartaer);
try
  [headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
catch me
  me
  %fprintf(1,'oops : error running sarta clear, look at error log %s \n',ugh);
  %keyboard
  disp('in get_sarta_clear.m')
  which mktemp
  which get_sarta_clear
  fprintf(1,'oops : error running sarta clear, look at ip/op rtp files %s %s \n',fip,fop)
  fprintf(1,'oops : error running sarta clear, look at error log - here is a snapshot %s \n',ugh2);
  morer = ['!more ' ugh2];
  eval(morer)    
  error('woof! try again!')
end

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh ' ' ugh2]; eval(rmer);

toc

hlayers = headRX2;
players = profRX2;

profOUT.sarta_rclearcalc = profRX2.rcalc;
