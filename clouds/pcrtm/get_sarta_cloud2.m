%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

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

if isfield(p,'co2ppm')
  p = rmfield(p,'co2ppm');
end
p.co2ppm = co2;    %% for cloud, use PCRTM values of co2 for debug purposes

p2junk = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

rtpwrite(fip,h,ha,p2junk,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
  eval(klayerser);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh];
  eval(sartaer);
[headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh]; eval(rmer);

[ppmvLAY2,ppmvAVG2,ppmvMAX2,pavgLAY2,tavgLAY2] = layers2ppmv(headRX,profRX,1:length(profRX.stemp),2);

figure(1);
plot(ppmvLAY2(40,:)); title('CO2 ppmv')

