%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

disp(' ')

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
ugh = mktemp('ugh');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(p,'co2ppm')
  p = rmfield(p,'co2ppm');
end

%% OLD 
% p.co2ppm = co2;    %% for cloud, use PCRTM values of co2 for debug purposes, ORIG
% p2junk  = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);   %% this has all the extra debugging cloud info
                                                           %% COMMENT OUT as we are adding in CO2 profile from XiuHong
							   
%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEW
if isfield(p,'co2ppm')
  p = rmfield(p,'co2ppm');
end

p_sarta_add_co2_ch4

p2junk  = driver_sarta_cloud_rtp(hxjunk,ha,pxjunk,pa,run_sarta);   %% this has all the extra debugging cloud info with NEW CO2 profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

profRX2 = p2junk;                                        %% new, since we HAVE run sarta cloud!!!

if run_sarta.clear > 0 & run_sarta.cloud > 0
  disp('ran off SARTA cloudy and SARTA clear .... in get_sarta_cloud2 (profRX2) and in get_sarta_clear2 (profRX)')
  fprintf(1,'check CLEAR calcs : difference sum(sum(p2junk.sarta_rclearcalc-profRX.rcalc)) = %8.6f \n',sum(sum(p2junk.sarta_rclearcalc-profRX.rcalc)))
end

%{
rtpwrite(fip,h,ha,p2junk,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
  eval(klayerser);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh];
  eval(sartaer);
try
  [headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
catch me
  me
  fprintf(1,'oops : error running sarta cloudy, look at error log %s \n',ugh2);
  error('woof! try again!')
end
%}

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh]; eval(rmer);

[ppmvLAY2,ppmvAVG2,ppmvMAX2,pavgLAY2,tavgLAY2] = layers2ppmv(headRX,profRX,1:length(profRX.stemp),2);

%figure(1); plot(ppmvLAY2(40,:)); title('CO2 ppmv')

