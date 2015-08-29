pcrtm_p = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/Atm_prof/lev-101_nMol-6/pbnd.dat');
pcrtm = load('/home/sergio/PCRTM_XIANGLEI/PCRTM_V2.1/code_changed/InputDir/par_constant.dat');
%% should be CO2 N2O CO CH4

hxjunk = h;
pxjunk = p;

[mmjunk,nnjunk] = size(pxjunk.plevs);
fprintf(1,'    >> size of plevs before adding in co2/ch4 = %5i x %5i \n',mmjunk,nnjunk)

hxjunk.ngas  = 6;
hxjunk.glist = [ 1  2  3  4  5  6]';
hxjunk.gunit = [21 10 21 10 10 10]';  %% mr g/g = 21,  ppmv = 10

for iijunk = 1 : length(pxjunk.stemp)
  nlevs = pxjunk.nlevs(iijunk);
  boo   = pxjunk.plevs(:,iijunk);

  %% co2 = 385.84 ppm at GND
  %junk = interp1(log10(pcrtm_p),pcrtm(:,1),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap') * co2(iijunk)/385.14;
  junk = interp1(log10(pcrtm_p),pcrtm(:,1),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap') * sarta_gas_2_6.co2/385.848;  
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_2(:,iijunk) = junky;

  %% n2o
  junk = interp1(log10(pcrtm_p),pcrtm(:,2),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_4(:,iijunk) = junky;

  %% co
  junk = interp1(log10(pcrtm_p),pcrtm(:,3),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_5(:,iijunk) = junky;

  %% ch4 = 1.843 ppm at GND
  %junk = interp1(log10(pcrtm_p),pcrtm(:,4),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap');
  junk = interp1(log10(pcrtm_p),pcrtm(:,4),log10(pxjunk.plevs(1:nlevs,iijunk)),'linear','extrap') * sarta_gas_2_6.ch4/1.843;  
  junky = ones(length(boo),1) * -9999;
  junky(1:nlevs) = junk;
  pxjunk.gas_6(:,iijunk) = junky;
  
end
