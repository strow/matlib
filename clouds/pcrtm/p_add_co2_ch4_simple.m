pcrtm_p = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/Atm_prof/lev-101_nMol-6/pbnd.dat');
pcrtm = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/par_constant.dat');
%% should be CO2 N2O CO CH4

hxjunk = h;
pxjunk = p;

hxjunk.ngas = 4;
hxjunk.glist = [ 1  2  3  6]';
hxjunk.gunit = [21 10 21 10]';  %% mr g/g = 21,  ppmv = 10

for iijunk = 1 : length(pxjunk.stemp)
  nlevs = pxjunk.nlevs(iijunk);
  boo   = pxjunk.plevs(:,iijunk);

  %% co2 = 385.84 ppm at GND
  %junk = interp1(log10(pcrtm_p),pcrtm(:,1),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap') * co2(iijunk)/385.14;
  junk = interp1(log10(pcrtm_p),pcrtm(:,1),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap') * sarta_gas_2_6.co2/385.14;  
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_2(:,iijunk) = junky;

  %% ch4 = 1.843 ppm at GND
  %junk = interp1(log10(pcrtm_p),pcrtm(:,4),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junk = interp1(log10(pcrtm_p),pcrtm(:,4),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap') * sarta_gas_2_6.ch4/1.843;  
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_6(:,iijunk) = junky;
  
end
