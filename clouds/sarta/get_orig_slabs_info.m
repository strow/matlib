function [orig_slabs,p] = get_orig_slabs_info(p,run_sarta);

orig_slabs = [];

iFound = 0;
if isfield(p,'cfrac')
  orig_slabs.cfrac = p.cfrac;
  iFound = iFound + 1;
end
if isfield(p,'ctype')
  orig_slabs.ctype = p.ctype;
  iFound = iFound + 1;  
end
if isfield(p,'cprtop')
  orig_slabs.cprtop = p.cprtop;
  iFound = iFound + 1;  
end
if isfield(p,'cprbot')
  orig_slabs.cprbot = p.cprbot;
  iFound = iFound + 1;  
end
if isfield(p,'cngwat')
  orig_slabs.cngwat = p.cngwat;
  iFound = iFound + 1;  
end
if isfield(p,'cpsize')
  orig_slabs.cpsize = p.cpsize;
  iFound = iFound + 1;  
end

if isfield(p,'cfrac2')
  orig_slabs.cfrac2 = p.cfrac2;
  iFound = iFound + 1;  
end
if isfield(p,'ctype2')
  orig_slabs.ctype2 = p.ctype2;
  iFound = iFound + 1;  
end
if isfield(p,'cprtop2')
  orig_slabs.cprtop2 = p.cprtop2;
  iFound = iFound + 1;  
end
if isfield(p,'cprbot2')
  orig_slabs.cprbot2 = p.cprbot2;
  iFound = iFound + 1;  
end
if isfield(p,'cngwat2')
  orig_slabs.cngwat2 = p.cngwat2;
  iFound = iFound + 1;  
end
if isfield(p,'cpsize2')
  orig_slabs.cpsize2 = p.cpsize2;
  iFound = iFound + 1;  
end

if isfield(p,'cfrac12')
  orig_slabs.cfrac12 = p.cfrac12;
  iFound = iFound + 1;  
end

if run_sarta.talk == 1
  if iFound > 0
    fprintf(1,'   found %2i of total 13 possible cloud fields in orig structure \n',iFound)
    fprintf(1,'   saving them so you can compare old vs new cloud slabs \n\n')
  elseif iFound == 0 & ~isfield(p,'tcc')
    disp('did not find any cloud fslb ields, including cfrac OR tcc !!!! in orig input structure, returning "orig_slabs = []" ...');
    disp('this means need to initialize cfracs ... and get pretty bad biases overall ... ugh')
    disp(' ')
  elseif iFound == 0 & isfield(p,'tcc')
    disp('did not find any cloud slab fields, but found tcc !!!! in orig input structure, returning "orig_slabs = []" ...');
    disp('this means need to initialize cfracs ... ')
    disp(' ')
  end
end

if ~isfield(p,'cfrac') & ~isfield(p,'tcc')
  %% need random cfracs
  if run_sarta.talk == 1
    disp('>>>>>>>> warning : need random cfracs, and NO tcc .... initializing')
  end
  %% want to make sure there are NO zeros cfrac
  p.cfrac = 0.50*(rand(size(p.stemp)) + rand(size(p.stemp)));
elseif ~isfield(p,'cfrac') & isfield(p,'tcc')
  %% need random cfracs
  if run_sarta.talk == 1
    disp('>>>>>>>> warning : need random cfracs, and YES tcc .... initializing')
  end
  %% want to make sure there are NO zeros cfrac
  p.cfrac = p.tcc;
end

if isfield(p,'tcc') & run_sarta.tcc > 0
  if run_sarta.talk == 1  
    disp(' woohoo : you have field "tcc" so resetting "cfrac" with this!!!')
  end
  p.cfrac = p.tcc;
elseif isfield(p,'tcc') & run_sarta.tcc < 0
  if run_sarta.talk == 1
    disp(' unfortunately run_sarta.tcc < 0 while you DO have p.tcc')
  end
elseif ~isfield(p,'tcc') & run_sarta.tcc > 0
  if run_sarta.talk == 1
    disp(' unfortunately run_sarta.tcc > 0 but you do NOT have p.tcc')
  end
  error('code now requires input field tcc')
elseif ~isfield(p,'tcc')
  error('code now requires input field tcc')  
end

disp(' ')
