foo = zeros(size(p.stemp));
boo = find(p.ctype == 101);  
  if length(boo) > 0
    foo(boo) = p.cngwat(boo);
  end
boo = find(p.ctype2 == 101);
  if length(boo) > 0
    foo(boo) = p.cngwat2(boo);
  end
wateramtS = [wateramtS foo(ix)];
wateramtP = [wateramtP p.pcrtm_waterOD(ix)];

foo = zeros(size(p.stemp));
boo = find(p.ctype == 201);  
  if length(boo) > 0
    foo(boo) = p.cngwat(boo);
  end
boo = find(p.ctype2 == 201);
  if length(boo) > 0
    foo(boo) = p.cngwat2(boo);
  end
iceamtS = [iceamtS foo(ix)];
iceamtP = [iceamtP p.pcrtm_iceOD(ix)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foo = zeros(size(p.stemp));
boo = find(p.ctype == 101);  
  if length(boo) > 0
    foo(boo) = p.cpsize(boo);
  end
boo = find(p.ctype2 == 101);
  if length(boo) > 0
    foo(boo) = p.cpsize2(boo);
  end
waterszeS = [waterszeS foo(ix)];
waterszeP = [waterszeP p.pcrtm_waterDME(ix)];

foo = zeros(size(p.stemp));
boo = find(p.ctype == 201);  
  if length(boo) > 0
    foo(boo) = p.cpsize(boo);
  end
boo = find(p.ctype2 == 201);
  if length(boo) > 0
    foo(boo) = p.cpsize2(boo);
  end
iceszeS = [iceszeS foo(ix)];
iceszeP = [iceszeP p.pcrtm_iceDME(ix)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foo = zeros(size(p.stemp));
boo = find(p.ctype == 101);  
  if length(boo) > 0
    foo(boo) = p.cprtop(boo);
  end
boo = find(p.ctype2 == 101);
  if length(boo) > 0
    foo(boo) = p.cprtop2(boo);
  end
watertopS = [watertopS foo(ix)];
watertopP = [watertopP p.pcrtm_waterCTOP(ix)];

foo = zeros(size(p.stemp));
boo = find(p.ctype == 201);  
  if length(boo) > 0
    foo(boo) = p.cprtop(boo);
  end
boo = find(p.ctype2 == 201);
  if length(boo) > 0
    foo(boo) = p.cprtop2(boo);
  end
icetopS = [icetopS foo(ix)];
icetopP = [icetopP p.pcrtm_iceCTOP(ix)];
