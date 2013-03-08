function [p1,oo] = only_waterORice_cloud(h,p0,waterORice)

%% p0 is input profile
%% p1 is output profile, keeping only "oo" from the input profile
%% depending on "waterORice"

if waterORice > 0
  disp(' >>>>>>>>>>>>> WARNING : turning off ice clouds <<<<<<<<<<<<<<<<')
  %% water is 101, ice is 201
  %% keep only water slab clouds %%
  oo1 = find(p0.ctype  == 201 & p0.ctype2 == 101); %% so 2 = water, 1 = ice
  oo2 = find(p0.ctype2 == 201 & p0.ctype  == 101); %% so 2 = ice,   1 = water
  %% these are the ones we want yay yay yay
  oo = union(oo1,oo2);
  [hx,p1x] = subset_rtp_allcloudfields(h,p0,[],[],oo);

  p1 = p1x;
  p1.ciwc = 0 * p1.ciwc;        %% set ice to 0
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumwater(ii) = trapz(p1.plevs(1:nlevs,ii),p1.clwc(1:nlevs,ii));
  end

  oo1 = find(p1.ctype  == 201 & p1.ctype2 == 101); %% so 2 = water, 1 = ice
  oo2 = find(p1.ctype2 == 201 & p1.ctype  == 101); %% so 2 = ice,   1 = water

  for ii = 1 : length(oo1)
    p1.cfrac(oo1(ii))   = 0;
    p1.cfrac12(oo1(ii)) = 0;
    p1.cngwat(oo1(ii))  = 0;
    p1.ctype(oo1(ii))   = -9999;
    cT = p1.cprtop2(oo1(ii));
    cB = p1.cprbot2(oo1(ii));
    p1.clwc(:,oo1(ii)) = 0;
    xx = find(p1.plevs(:,oo1(ii)) >= cT & p1.plevs(:,oo1(ii)) <= cB);
    p1.clwc(xx,oo1(ii)) = sumwater(oo1(ii))/(cB-cT);
  end
  for ii = 1 : length(oo2)
    p1.cfrac2(oo2(ii))  = 0;
    p1.cfrac12(oo2(ii)) = 0;
    p1.cngwat2(oo2(ii)) = 0;
    p1.ctype2(oo2(ii))  = -9999;
    cT = p1.cprtop(oo2(ii));
    cB = p1.cprbot(oo2(ii));
    p1.clwc(:,oo2(ii)) = 0;
    xx = find(p1.plevs(:,oo2(ii)) >= cT & p1.plevs(:,oo2(ii)) <= cB);
    p1.clwc(xx,oo2(ii)) = sumwater(oo2(ii))/(cB-cT);
  end
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumwater2(ii) = trapz(p1.plevs(1:nlevs,ii),p1.clwc(1:nlevs,ii));
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  for ii = 1 : length(oo1)
    p1.cfrac(oo1(ii))   = p1.cfrac2(oo1(ii));
    p1.cngwat(oo1(ii))  = p1.cngwat2(oo1(ii));
    p1.cprtop(oo1(ii))  = p1.cprtop2(oo1(ii));
    p1.cprbot(oo1(ii))  = p1.cprbot2(oo1(ii));
    p1.cpsize(oo1(ii))  = p1.cpsize2(oo1(ii));
    p1.ctype(oo1(ii))   = 101;
    p1.cfrac12(oo1(ii)) = 0;
    p1.cfrac2(oo1(ii))  = 0;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waterORice < 0
  disp(' >>>>>>>>>>>>> WARNING : turning off water clouds <<<<<<<<<<<<<<<<')
  %% water is 101, ice is 201
  %% keep only ice slab clouds %%
  oo1 = find(p0.ctype  == 101 & p0.ctype2 == 201); %% so 1 = water, 2 = ice
  oo2 = find(p0.ctype2 == 101 & p0.ctype  == 201); %% so 1 = ice,   2 = water
  %% these are the ones we want yay yay yay
  oo = union(oo1,oo2);
  [hx,p1x] = subset_rtp_allcloudfields(h,p0,[],[],oo);

  p1 = p1x;
  p1.clwc = 0 * p1.clwc;        %% set water to 0
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumice(ii) = trapz(p1.plevs(1:nlevs,ii),p1.ciwc(1:nlevs,ii));
  end

  oo1 = find(p1.ctype  == 101 & p1.ctype2 == 201); %% so 1 = water, 2 = ice
  oo2 = find(p1.ctype2 == 101 & p1.ctype  == 201); %% so 1 = ice,   2 = water

  for ii = 1 : length(oo1)
    p1.cfrac(oo1(ii))   = 0;
    p1.cfrac12(oo1(ii)) = 0;
    p1.cngwat(oo1(ii))  = 0;
    p1.ctype(oo1(ii))   = -9999;
    cT = p1.cprtop2(oo1(ii));
    cB = p1.cprbot2(oo1(ii));
    p1.ciwc(:,oo1(ii)) = 0;
    xx = find(p1.plevs(:,oo1(ii)) >= cT & p1.plevs(:,oo1(ii)) <= cB);
    p1.ciwc(xx,oo1(ii)) = sumice(oo1(ii))/(cB-cT);
  end
  for ii = 1 : length(oo2)
    p1.cfrac2(oo2(ii))  = 0;
    p1.cfrac12(oo2(ii)) = 0;
    p1.cngwat2(oo2(ii)) = 0;
    p1.ctype2(oo2(ii))  = -9999;
    cT = p1.cprtop(oo2(ii));
    cB = p1.cprbot(oo2(ii));
    p1.ciwc(:,oo2(ii)) = 0;
    xx = find(p1.plevs(:,oo2(ii)) >= cT & p1.plevs(:,oo2(ii)) <= cB);
    p1.ciwc(xx,oo2(ii)) = sumice(oo2(ii))/(cB-cT);
  end
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumice2(ii) = trapz(p1.plevs(1:nlevs,ii),p1.ciwc(1:nlevs,ii));
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  for ii = 1 : length(oo1)
    p1.cfrac(oo1(ii))   = p1.cfrac2(oo1(ii));
    p1.cngwat(oo1(ii))  = p1.cngwat2(oo1(ii));
    p1.cprtop(oo1(ii))  = p1.cprtop2(oo1(ii));
    p1.cprbot(oo1(ii))  = p1.cprbot2(oo1(ii));
    p1.cpsize(oo1(ii))  = p1.cpsize2(oo1(ii));
    p1.ctype(oo1(ii))   = 201;
    p1.cfrac12(oo1(ii)) = 0;
    p1.cfrac2(oo1(ii))  = 0;

  end

end
