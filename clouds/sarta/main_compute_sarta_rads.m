%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  adds CO2 profile, if needed removes ice or water clds, runs sarta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iCompareSlabs = -1;
iCompareSlabs = +1;
if iCompareSlabs > 0
  %% compare pINPUT and p
  dn = 0 : 0.1 : 1;
  figure(1); clf;
    nx1 = hist(pINPUT.cfrac,dn);  ny1= hist(prof.cfrac,dn);
    nx2 = hist(pINPUT.cfrac2,dn); ny2= hist(prof.cfrac2,dn);
    nx12 = hist(pINPUT.cfrac12,dn); ny12= hist(prof.cfrac12,dn);      
    plot(dn,nx1,'bo-',dn,nx2,'bx--',dn,ny1,'r',dn,ny2,'r--',dn,nx12,'ko-',dn,ny12,'k--','linewidth',2)
    hl = legend('pIN cfrac1','pIN cfrac2','prof cfrac1','prof cfrac2','pIN cfrac12','prof cfrac12');
    set(hl,'fontsize',10);

  dn = 0 : 50 : 1000;
  figure(2); clf;
    nx1 = hist(pINPUT.cprtop,dn);  ny1= hist(prof.cprtop,dn);
    nx2 = hist(pINPUT.cprtop2,dn); ny2= hist(prof.cprtop2,dn);  
    plot(dn,nx1,'bo-',dn,nx2,'bx--',dn,ny1,'r',dn,ny2,'r--','linewidth',2)
    hl = legend('pIN cprtop1','pIN cprtop2','prof cprtop1','prof cprtop2');
    set(hl,'fontsize',10);

  dn = 0 : 50 : 1000;
  figure(3); clf;
    nx1 = hist(pINPUT.cprbot,dn);  ny1= hist(prof.cprbot,dn);
    nx2 = hist(pINPUT.cprbot2,dn); ny2= hist(prof.cprbot2,dn);  
    plot(dn,nx1,'bo-',dn,nx2,'bx--',dn,ny1,'r',dn,ny2,'r--','linewidth',2)
    hl = legend('pIN cprbot1','pIN cprbot2','prof cprbot1','prof cprbot2');
    set(hl,'fontsize',10);

  dn = 0 : 2 : 500;
  figure(4); clf;
    nx1 = hist(pINPUT.cngwat,dn);  ny1= hist(prof.cngwat,dn);
    nx2 = hist(pINPUT.cngwat2,dn); ny2= hist(prof.cngwat2,dn);  
    loglog(dn,nx1,'bo-',dn,nx2,'bx--',dn,ny1,'r',dn,ny2,'r--','linewidth',2)
    hl = legend('pIN cngwat1','pIN cngwat2','prof cngwat1','prof cngwat2');
    set(hl,'fontsize',10);
    ax = axis; axis([2 500 ax(3) ax(4)])
    
  dn = 0 : 5 : 200;
  figure(5); clf;
    nx1 = hist(pINPUT.cpsize,dn);  ny1= hist(prof.cpsize,dn);
    nx2 = hist(pINPUT.cpsize2,dn); ny2= hist(prof.cpsize2,dn);  
    semilogy(dn,nx1,'bo-',dn,nx2,'bx--',dn,ny1,'r',dn,ny2,'r--','linewidth',2)
    hl = legend('pIN cpsize1','pIN cpsize2','prof cpsize1','prof cpsize2');
    set(hl,'fontsize',10);

  junk = [mean(pINPUT.cfrac-prof.cfrac) mean(pINPUT.cfrac2-prof.cfrac2) mean(pINPUT.cfrac12-prof.cfrac12)];
  fprintf(1,'checking profile : cfrac1,cfrac2,cfrac12 diffs = %8.6f %8.6f %8.6f \n',junk);

  oo11 = find(pINPUT.ctype  == 101);   oo12 = find(pINPUT.ctype  == 201);
  oo21 = find(pINPUT.ctype2 == 101);   oo22 = find(pINPUT.ctype2 == 201);
  oo11 = length(oo11); oo12 = length(oo12); oo21 = length(oo21); oo22 = length(oo22);
  junk = [oo11 oo12 oo21 oo22];
  fprintf(1,'pINPUT : number of ctype=101/201 ctype2=101/201 = %4i %4i %4i %4i \n',junk);

  oo11 = find(prof.ctype  == 101);   oo12 = find(prof.ctype  == 201);
  oo21 = find(prof.ctype2 == 101);   oo22 = find(prof.ctype2 == 201);
  oo11 = length(oo11); oo12 = length(oo12); oo21 = length(oo21); oo22 = length(oo22);  
  junk = [oo11 oo12 oo21 oo22];
  fprintf(1,'prof   : number of ctype=101/201 ctype2=101/201 = %4i %4i %4i %4i \n',junk);  

  fprintf(1,'checking profile : stemp diff = %8.6f \n',nansum(pINPUT.stemp-prof.stemp))
  fprintf(1,'checking profile : plevs diff = %8.6f \n',nansum(nansum(pINPUT.plevs-prof.plevs)))  
  fprintf(1,'checking profile : ptemp diff = %8.6f \n',nansum(nansum(pINPUT.ptemp-prof.ptemp)))  
  fprintf(1,'checking profile : gas_1 diff = %8.6f \n',nansum(nansum(pINPUT.gas_1-prof.gas_1)))
  fprintf(1,'checking profile : gas_3 diff = %8.6f \n',nansum(nansum(pINPUT.gas_3-prof.gas_3)))    
  pause(0.1);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add on co2
prof_add_co2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now if wanted can turn off ice or water clouds!!!!!!!!
%% this is kinda doing what "driver_sarta_cloud_rtp_onecldtest.m" is meant to do
if run_sarta.waterORice ~= 0
  disp('run_sarta.waterORice ~= 0 ---->>> going to remove ice or water cld')
  [prof,index_kept] = only_waterORice_cloud(h,prof,run_sarta.waterORice);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  tic
  get_sarta_clear;
  toc
  prof.sarta_rclearcalc = profRX2.rcalc;
else
  disp('you did not ask for SARTA clear to be run; not changing prof.sarta_rclearcalc')  
end

if run_sarta.cloud > 0 
  disp('running SARTA cloud')
  tic
  get_sarta_cloud;
  toc
  prof.rcalc = profRX2.rcalc;
else
  disp('you did not ask for SARTA cloudy to be run; not changing prof.rcalc')
end
