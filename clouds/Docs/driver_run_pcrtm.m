% local running to test
% clustcmd -L driver_run_pcrtm.m 1
%
% otherwise when happy
% clustcmd -q medium -n 32  driver_run_pcrtm.m 1:24
%
% or
% clustcmd -q long_contrib -n 24  driver_run_pcrtm.m 1:24
% clustcmd -q medium -P test  -n 24  driver_run_pcrtm.m 1:24

%% make sure you have " pcrtm.in.tmp " in this rundir

addpath /asl/matlib/gribtools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil/
addpath /asl/matlib/science/
addpath /asl/matlib/h4tools/

addpath /strowdata1/shared/sergio/MATLABCODE/

addpath /home/sergio/MATLABCODE/matlib/clouds/pcrtm
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

fairs = instr_chans;

dd = 1;
dd = 2;
dd = 3;
dd = 4;
dd = 5;
dd = 6;
dd = 7;
dd = 8;
dd = 9;
dd = 10;

ddStart = 01; ddEnd = 05;
ddStart = 06; ddEnd = 10;
ddStart = 11; ddEnd = 15;
ddStart = 16; ddEnd = 20;
ddStart = 21; ddEnd = 25;
ddStart = 26; ddEnd = 31;

ddStart = 01; ddEnd = 01;
ddStart = 01; ddEnd = 10;
ddStart = 11; ddEnd = 20;
ddStart = 21; ddEnd = 31;

ddStart = 01; ddEnd = 10;

for dd = ddStart : ddEnd
  clear h ha p pa

  fmain = ['/asl/data/rtprod_airs/2012/05/' num2str(dd,'%02d') '/cld_ecm_41ch.airs_ctr.2012.05.' num2str(dd,'%02d') '.'];
  fin = [fmain  num2str(JOB,'%02d') '.rtp'];                                  %% oops, forget about Scotts new random cld frac
  fout = [fmain  num2str(JOB,'%02d') '_2378chansNEW.rtp'];                    %% fixed the random cld fracn when calling sarta
  fout = [fmain  num2str(JOB,'%02d') '_2378chansNEW_ncol0_1.rtp'];            %% fixed the random cld fracn when calling sarta
  fout = [fmain  num2str(JOB,'%02d') '_2378chansNEW_ncol0_1_csum_0p3.rtp'];   %% fixed the random cld fracn when calling sarta
  fout = [fmain  num2str(JOB,'%02d') '_2378chansNEW2_ncol0_1.rtp'];           %% oops, found a bug in ecmwf2sarta, testing

  ee = exist(fin);
  if ee > 0
    [h,ha,p,pa] = rtpread(fin);
    [h,ha,p,pa] = rtpgrow(h,ha,p,pa);

    robs1 = p.robs1;
    rcalc = p.rcalc;
    calflag = p.calflag;
    ichan = h.ichan;

    h.nchan = 2378;
    h.ichan = (1:2378)';
    h.vchan = fairs;

    p.robs1 = zeros(2378,length(p.stemp));
    p.robs1(ichan,:) = robs1;
    p.rcalc = zeros(2378,length(p.stemp));
    p.rcalc(ichan,:) = rcalc;
    p.calflag = zeros(2378,length(p.stemp));
    p.calflag(ichan,:) = calflag;

    run_sarta.clear = +1;
    run_sarta.cloud = +1;
    run_sarta.ncol0 = +1;    %%% instead of 50
%    run_sarta.cumsum = 0.3;  %%% this is for the sarta run

    tic
    p1 = driver_pcrtm_cloud_rtp(h,ha,p,pa,run_sarta);
    toc

    rtpwrite(fout,h,ha,p1,pa);

  end  %% if ee > 0
end    %% for dd loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
