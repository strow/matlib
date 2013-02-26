addpath /asl/matlab/aslutil/

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

fip = mktemp('temp.ip.rtp');
fop = mktemp('temp.op.rtp');
frp = mktemp('temp.rp.rtp');
ugh = mktemp('ugh');

if isfield(p,'co2ppm')
  p = rmfield(p,'co2ppm');
end
p.co2ppm = co2;    %% for clear, use PCRTM values of co2 for debug purposes

oldrtpwrite(fip,h,ha,p,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
  eval(klayerser);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh];
  eval(sartaer);
[headRX hattrR profRX pattrR] = oldrtpread(frp);
rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh]; eval(rmer);

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY] = layers2ppmv(headRX,profRX,1:length(profRX.stemp),2);

figure(1);
plot(ppmvLAY(40,:)); title('CO2 ppmv')

