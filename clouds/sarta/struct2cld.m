function   [xcprtopS1,xcngwatS1,xcfracS1,xcpsize1,xsumODS1,xcprtopS2,xcngwatS2,xcfracS2,xcpsize2,xsumODS2] = ...
      struct2cld(xpS,xcprtopS1,xcngwatS1,xcfracS1,xcpsize1,xsumODS1,xcprtopS2,xcngwatS2,xcfracS2,xcpsize2,xsumODS2);

boo1 = ones(size(xpS.rtime))*nan; 
boo2 = ones(size(xpS.rtime))*nan; 
boo3 = ones(size(xpS.rtime))*nan; 
boo4 = ones(size(xpS.rtime))*nan; 
boo5 = ones(size(xpS.rtime))*nan;

oo = find(xpS.ctype == 201);  %% ice 
  boo1(oo) = xpS.cprtop(oo);   
  boo2(oo) = xpS.cngwat(oo);   
  boo3(oo) = xpS.cfrac(oo);    
  boo4(oo) = xpS.cpsize(oo);                  

  if isfield(xpS,'sarta_lvlODice')
    boo5(oo) = nansum(xpS.sarta_lvlODice(:,oo));
  end

oo = find(xpS.ctype2 == 201);  %% ice 
  boo1(oo) = xpS.cprtop2(oo); 
    xcprtopS2 = [xcprtopS2 boo1];   
  boo2(oo) = xpS.cngwat2(oo); 
    xcngwatS2 = [xcngwatS2 boo2];   
  boo3(oo) = xpS.cfrac2(oo);  
    xcfracS2  = [xcfracS2 boo3];
  boo4(oo) = xpS.cpsize2(oo); 
    xcpsize2 = [xcpsize2 boo5];
  if isfield(xpS,'sarta_lvlODice')
    boo5(oo) = boo5(oo) + nansum(xpS.sarta_lvlODice(:,oo)); 
    xsumODS2 = [xsumODS2 boo5];
  else
    xsumODS2 = [xsumODS2 boo5];
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boo1 = ones(size(xpS.rtime))*nan; 
boo2 = ones(size(xpS.rtime))*nan; 
boo3 = ones(size(xpS.rtime))*nan; 
boo4 = ones(size(xpS.rtime))*nan; 
boo5 = ones(size(xpS.rtime))*nan;

oo = find(xpS.ctype == 101);  %% water
  boo1(oo) = xpS.cprtop(oo);   
  boo2(oo) = xpS.cngwat(oo);   
  boo3(oo) = xpS.cfrac(oo);    
  boo4(oo) = xpS.cpsize(oo);                  
  if isfield(xpS,'sarta_lvlODwater')
    boo5(oo) = nansum(xpS.sarta_lvlODwater(:,oo));
  end

oo = find(xpS.ctype2 == 101);  %% water 
  boo1(oo) = xpS.cprtop2(oo); 
    xcprtopS1 = [xcprtopS1 boo1];   
  boo2(oo) = xpS.cngwat2(oo); 
    xcngwatS1 = [xcngwatS1 boo2];   
  boo3(oo) = xpS.cfrac2(oo);  
    xcfracS1  = [xcfracS1 boo3];
  boo4(oo) = xpS.cpsize2(oo); 
    xcpsize1 = [xcpsize1 boo5];
  if isfield(xpS,'sarta_lvlODwater')
    boo5(oo) = boo5(oo) + nansum(xpS.sarta_lvlODwater(:,oo)); 
    xsumODS1 = [xsumODS1 boo5];
  else
    xsumODS1 = [xsumODS1 boo5];
  end
