function rcalcdiff = compare_profstructs_allprofs(h,p00,p10);

%% simple routine to compare profilesinput rtp structs  p1 and p2

plot(sum(p00.clwc),p00.cngwat,'.',sum(p10.clwc),p10.cngwat,'ro',sum(p00.clwc),p00.cngwat2,'.',sum(p10.clwc),p10.cngwat2,'ro')
  title('cngwat vs WATER sum(clwc), cngwat2 vs WATER sum(clwc)'); 
  disp('ret to continue'); pause;

plot(sum(p00.ciwc),p00.cngwat,'.',sum(p10.ciwc),p10.cngwat,'ro',sum(p00.ciwc),p00.cngwat2,'.',sum(p10.ciwc),p10.cngwat2,'ro')
  title('cngwat vs ICE sum(ciwc), cngwat2 vs ICE sum(ciwc)'); 
  disp('ret to continue'); pause;

plot(sum(p00.clwc)+sum(p00.ciwc),p00.cngwat+p00.cngwat2,'.',sum(p10.clwc)+sum(p10.ciwc),p10.cngwat+p10.cngwat2,'ro')
  title('cngwat+cngwat2 vs sum(ciwc)+sum(clwc)'); 
  disp('ret to continue'); pause;

plot(sum(p00.clwc)-sum(p10.clwc),p10.cngwat,'ro')                
  title('cngwat vs WATER sum(clwc1-clwc2)'); 
  disp('ret to continue'); pause;

plot(sum(p00.clwc),p00.cngwat-p10.cngwat,'ro',sum(p00.clwc),p00.cngwat2-p10.cngwat2,'ro')
  title('cngwat diff and cngwat2 diff vs WATER sum(clwc1)'); 
  disp('ret to continue'); pause;

plot(p10.ciwc./p00.ciwc)  
  title('ICE ciwc ratio'); 
  disp('ret to continue'); pause;

plot(p10.clwc./p00.clwc)
  title('WATER clwc ratio'); 
  disp('ret to continue'); pause;

plot(1:length(p10.stemp),p10.cprtop-p00.cprtop,1:length(p10.stemp),p10.cprtop2-p00.cprtop2,'r')
  title('cprtop diff'); 
  disp('ret to continue'); pause;

plot(1:length(p10.stemp),p10.cprbot-p00.cprbot,1:length(p10.stemp),p10.cprbot2-p00.cprbot2,'r')
  title('cprbot diff'); 
  disp('ret to continue'); pause;

plot(1:length(p10.stemp),p10.cngwat-p00.cngwat,1:length(p10.stemp),p10.cngwat2-p00.cngwat2,'r')
  title('cngwat diff'); 
  disp('ret to continue'); pause;

plot(1:length(p10.stemp),p10.cpsize-p00.cpsize,1:length(p10.stemp),p10.cpsize2-p00.cpsize2,'r')
  title('cpsize diff'); 
  disp('ret to continue'); pause;

plot(1:length(p10.stemp),p10.cfrac-p00.cfrac,1:length(p10.stemp),p10.cfrac2-p00.cfrac2,'r',1:length(p10.stemp),p10.cfrac12-p00.cfrac12,'k')
  title('cfrac diff'); 
  disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(p00.stemp-p10.stemp)
  title('stemp diff'); 
  disp('ret to continue'); pause;

plot(p00.ptemp-p10.ptemp)
  title('ptemp diff'); 
  disp('ret to continue'); pause;

plot(p00.gas_1./p10.gas_1)
  title('gas1 ratio'); 
  disp('ret to continue'); pause;
