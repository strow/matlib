  h2 = head;
  h2.ngas = 4; 
  h2.glist = [  1  3 201 202]'; 
  h2.gunit = [21  21  21  21]'; 
  h2.gunit = [21  21  26  26]';  %%this is for "wet" mixing ratios 

  p2.gas_1    = prof.gas_1;  % water
  p2.gas_3    = prof.gas_3;  % water
  %p2.landfrac = prof.landfrac;
  %p2.salti    = prof.salti;
  p2.nlevs    = prof.nlevs;
  p2.plat     = profX.plat';
  p2.plon     = profX.plon';
  p2.stemp    = profX.stemp;
  p2.spres    = profX.spres;
  p2.plevs    = profX.plevs;
  p2.ptemp    = profX.ptemp;
  p2.gas_201  = profX.clwc;    % water cloud 
  p2.gas_202  = profX.ciwc;    % ice cloud 

  randstr = num2str(rand(1,1)*10000000);
  fip_cld2 = ['fcld.ip.rtp.' randstr];
  fop_cld2 = ['fcld.op.rtp.' randstr];
  ugh2     = ['ugh'          randstr];

  rtpwrite(fip_cld2,h2,[],p2,[]); 

  klayers = '/asl/packages/klayersV205/Bin/klayers_airs_trace '; 
  klayers = '/asl/packages/klayersV205/Bin/klayers_airs_trace_testme '; 
  klayers = '/asl/packages/klayersV205/Bin/klayers_airs_v5_testme '; 
  klayerser = ['! ' klayers ' nwant=6 listg=1,2,3,4,201,202 ']; 
  klayerser = [klayerser ' fin=' fip_cld2 ' fout=' fop_cld2 ' > ' ugh2 ]; 
  eval(klayerser); 

  [h3,ha3,p3,pa3] = rtpread(fop_cld2);

  rmer = ['!/bin/rm ' fip_cld2 ' ' fop_cld2 ' ' ugh2];
  eval(rmer);

  ice_from_p3         = sum(p3.gas_202);
  water_from_p3       = sum(p3.gas_201);
  slabice_from_prof   = (prof.cngwat);
  slabwater_from_prof = (prof.udef(11,:));
  clds = 1 : length(ice_from_p3);
  plot(clds,ice_from_p3,'b',  clds,slabice_from_prof,'b--',...
       clds,water_from_p3,'r',clds,slabwater_from_prof,'r--');
  title('blue/red = ice/water solid/dashed = scott/me');

