%% these are required but user needs to add them before using this code
%% addpath /asl/matlib/aslutil/

disp(' ')
fprintf(1,'  doing clear sky SARTA calcs ....\n')
fprintf(1,'    sarta clear = %s \n',run_sarta.sartaclear_code);  

klayers = run_sarta.klayers_code;
sarta   = run_sarta.sartaclear_code;

if ~exist(klayers,'file')
  error('klayers exec done not exist')
end
if ~exist(sarta,'file')
  error('sarta clear exec done not exist')
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
p.co2ppm = co2;           %% for clear, use PCRTM values of co2 for debug purposes, ORIG
rtpwrite(fip,h,ha,p,pa);  %% comment this out, as we are adding in CO2 profile from XIuhong

klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
eval(klayerser);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh];
eval(sartaer);
[xheadRX xhattrR xprofRX xpattrR] = rtpread(frp);  %% was [headRX hattrR profRX pattrR] = rtpread(frp);
[xppmvLAY,xppmvAVG,xppmvMAX,xpavgLAY,xtavgLAY] = layers2ppmv(xheadRX,xprofRX,1:length(xprofRX.stemp),2);
[yppmvLAY,yppmvAVG,yppmvMAX,ypavgLAY,ytavgLAY] = layers2ppmv(xheadRX,xprofRX,1:length(xprofRX.stemp),6);  
rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh]; eval(rmer);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEW
if isfield(p,'co2ppm')
  p = rmfield(p,'co2ppm');
end

[hxjunk,pxjunk] = p_sarta_add_co2_ch4(h,p,sarta_gas_2_6,p0ALL);

[mmjunk,nnjunk] = size(pxjunk.plevs);
fprintf(1,'    >> size of plevs after  adding in co2/ch4 = %5i x %5i \n',mmjunk,nnjunk)

rtpwrite(fip,hxjunk,ha,pxjunk,pa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
  eval(klayerser);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh];
  eval(sartaer);
[headRX hattrR profRX pattrR] = rtpread(frp);
rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh]; eval(rmer);

[ppmvLAY, ppmvAVG,  ppmvMAX, pavgLAY, tavgLAY] = layers2ppmv(headRX,profRX,1:length(profRX.stemp),2);
[ppmvLAY6, ppmvAVG6,ppmvMAX6,pavgLAY6,tavgLAY6] = layers2ppmv(headRX,profRX,1:length(profRX.stemp),6);

%{
keyboard
plot(xppmvLAY,xpavgLAY/100,'b',ppmvLAY,pavgLAY/100,'r');   axis([350 390 0 1000]); set(gca,'ydir','reverse')
plot(yppmvLAY,ypavgLAY/100,'b',ppmvLAY6,pavgLAY6/100,'r'); axis([1.5 2.2 0 1000]); set(gca,'ydir','reverse')
figure(1); plot(ppmvLAY(40,:)); title('CO2 ppmv')
%}

