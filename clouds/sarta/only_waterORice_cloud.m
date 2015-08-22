function [p1,index_cloud_selected,index_other] = only_waterORice_cloud(h,p0,waterORice)

%% p0 is input profile
%%      "waterORice" +1 find pixels with (one water slab & one ice slab), removes the ice slab cloud, rejigs water profile
%%                      find pixels with (one water slab)                                           , rejigs water profile
%%                   -1 find pixels with (one water slab & one ice slab), removes the water slab cloud, rejigs ice profile
%%                      find pixels with (one ice slab)                                             , rejigs ice profile
%%                    0 do   nothing
%% p1 is output profile,
%%   index_cloud_selected are those whose clouds were selected as the pixel had (one water and one ice slab)  only
%%   index_other          are those that are clear or both ice or both water
%% 
%% >>>>>>>>>>>
%% so basically if you want to test PCRTM vs SARTA SLAB cloud you should
%{
[h,ha,p,pa] = colocate_era_airs(yy,mm,dd,gg); or make_rtp_yada(yy,mm,dd,gg);
run_sarta.ncol0 = -1;      %%% instead of 50, use -1 for testing ONE cloud
run_sarta.waterORice = +1; %%% turns off ice, keeps water clouds
run_sarta.waterORice = -1; %%% turns off water, keeps ice clouds
[p1,index_cloud_selected,index_other] = only_waterORice_cloud(h,p0,run_sarta.waterORice); %% sub-select according to water or ice; re-jigs relevant cloud profile
p2 = driver_pcrtm_cloud_rtp(h,ha,p1,pa,run_sarta);                                        %% call sarta and pcrtm on the SLAB profiles !!!
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waterORice == 0
  disp(' >>>>>>>>>>>>> in only_waterORice_cloud.m, doing NOTHING (keeping ice/water clds) <<<<<<<<<<<<<<<<')
  p1 = p0;
  index_cloud_selected = [];
  index_other          = [];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(p0.ctype == 101) == 0 & length(p0.ctype == 201) == 0 & length(p0.ctype2 == 101) == 0 & length(p0.ctype2 == 201) == 0
  error('only_waterORice_cloud : assumes you have already turned profiles into slabs')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waterORice > 0
  disp(' >>>>>>>>>>>>> SARTA : turning off ice clouds <<<<<<<<<<<<<<<<')
  disp(' >>> keeping only water clouds ==> less profiles out than in <<<')  
  
  %% water is 101, ice is 201
  %% keep only water slab clouds %%
  index_cloud_selected1 = find(p0.ctype  == 201 & p0.ctype2 == 101); %% so 2 = water, 1 = ice
  index_cloud_selected2 = find(p0.ctype2 == 201 & p0.ctype  == 101); %% so 2 = ice,   1 = water
  index_cloud_selected3 = find(p0.ctype2 == 101 & p0.ctype   < 101); %% so 2 = water and 1 is not an ice cloud (ie clear)
  index_cloud_selected4 = find(p0.ctype  == 101 & p0.ctype2  < 101); %% so 2 is not an ice cloud (ie clear), 1 = water
  index_cloud_selected5 = find(p0.ctype  == 101 & p0.ctype2 == 101); %% so 1,2 are both water

  %% these are the ones we want yay yay yay
  index_cloud_selected = union(union(union(union(index_cloud_selected1,index_cloud_selected2),index_cloud_selected3),index_cloud_selected4),index_cloud_selected5);

  %% these are the ones we do not want bah bah bah
  index_other = 1 : length(p0.stemp);
  index_other = setdiff(index_other,index_cloud_selected); %% either clear or both water or both ice

  junk = [length(index_cloud_selected1) length(index_cloud_selected2) length(index_cloud_selected3) length(index_cloud_selected4) length(index_cloud_selected5) length(index_other)];
  fprintf(1,'number of I/W W/I X/W W/X W/W X/X = %4i %4i %4i %4i %4i %4i \n',junk)
  
  [hx,p1x] = subset_rtp_allcloudfields(h,p0,[],[],index_cloud_selected);
  if length(index_other) > 0
    [hx,p1z] = subset_rtp_allcloudfields(h,p0,[],[],index_other);
    p1z.ciwc = 0 * p1z.clwc;      %% set ice to 0
    p1z.clwc = 0 * p1z.clwc;      %% set water to 0, since it was not selected anyway!!! ???
  end
  
  p1 = p1x;
  p1.ciwc = 0 * p1.ciwc;        %% set ice to 0
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumwater(ii) = trapz(p1.plevs(1:nlevs,ii),p1.clwc(1:nlevs,ii));
  end

  index_cloud_selected1 = find(p1.ctype  == 201 & p1.ctype2 == 101); %% so 2 = water, 1 = ice
  index_cloud_selected2 = find(p1.ctype2 == 201 & p1.ctype  == 101); %% so 2 = ice,   1 = water

  for ii = 1 : length(index_cloud_selected1)
    p1.cfrac(index_cloud_selected1(ii))   = 0;
    p1.cfrac12(index_cloud_selected1(ii)) = 0;
    p1.cngwat(index_cloud_selected1(ii))  = 0;
    p1.ctype(index_cloud_selected1(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere    
    cT = p1.cprtop2(index_cloud_selected1(ii));
    cB = p1.cprbot2(index_cloud_selected1(ii));
    p1.clwc(:,index_cloud_selected1(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected1(ii)) >= cT & p1.plevs(:,index_cloud_selected1(ii)) <= cB);
    p1.clwc(xx,index_cloud_selected1(ii)) = sumwater(index_cloud_selected1(ii))/(cB-cT);
  end
  
  for ii = 1 : length(index_cloud_selected2)
    p1.cfrac2(index_cloud_selected2(ii))  = 0;
    p1.cfrac12(index_cloud_selected2(ii)) = 0;
    p1.cngwat2(index_cloud_selected2(ii)) = 0;
    p1.ctype2(index_cloud_selected2(ii))  = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(index_cloud_selected2(ii));
    cB = p1.cprbot(index_cloud_selected2(ii));
    p1.clwc(:,index_cloud_selected2(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected2(ii)) >= cT & p1.plevs(:,index_cloud_selected2(ii)) <= cB);
    p1.clwc(xx,index_cloud_selected2(ii)) = sumwater(index_cloud_selected2(ii))/(cB-cT);
    p1.cfrac(index_cloud_selected2(ii))  = 1;
  end

  index_cloud_selected3 = find(p1.ctype2 == 101 & p1.ctype  < 101); %% so 2 = water and 1 is not an ice cloud (ie clear)
  index_cloud_selected4 = find(p1.ctype  == 101 & p1.ctype2 < 101); %% so 2 is not an ice cloud (ie clear), 1 = water

  for ii = 1 : length(index_cloud_selected3)
    p1.cfrac(index_cloud_selected3(ii))   = 0;
    p1.cfrac12(index_cloud_selected3(ii)) = 0;
    p1.cngwat(index_cloud_selected3(ii))  = 0;
    p1.ctype(index_cloud_selected3(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere    
    cT = p1.cprtop2(index_cloud_selected3(ii));
    cB = p1.cprbot2(index_cloud_selected3(ii));
    p1.clwc(:,index_cloud_selected3(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected3(ii)) >= cT & p1.plevs(:,index_cloud_selected3(ii)) <= cB);
    p1.clwc(xx,index_cloud_selected3(ii)) = sumwater(index_cloud_selected3(ii))/(cB-cT);
  end
  
  for ii = 1 : length(index_cloud_selected4)
    p1.cfrac2(index_cloud_selected4(ii))  = 0;
    p1.cfrac12(index_cloud_selected4(ii)) = 0;
    p1.cngwat2(index_cloud_selected4(ii)) = 0;
    p1.ctype2(index_cloud_selected4(ii))  = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(index_cloud_selected4(ii));
    cB = p1.cprbot(index_cloud_selected4(ii));
    p1.clwc(:,index_cloud_selected4(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected4(ii)) >= cT & p1.plevs(:,index_cloud_selected4(ii)) <= cB);
    p1.clwc(xx,index_cloud_selected4(ii)) = sumwater(index_cloud_selected4(ii))/(cB-cT);
    p1.cfrac(index_cloud_selected4(ii))  = 1;
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  for ii = 1 : length(index_cloud_selected1)
    p1.cfrac(index_cloud_selected1(ii))   = p1.cfrac2(index_cloud_selected1(ii));
    p1.cngwat(index_cloud_selected1(ii))  = p1.cngwat2(index_cloud_selected1(ii));
    p1.cprtop(index_cloud_selected1(ii))  = p1.cprtop2(index_cloud_selected1(ii));
    p1.cprbot(index_cloud_selected1(ii))  = p1.cprbot2(index_cloud_selected1(ii));
    p1.cpsize(index_cloud_selected1(ii))  = p1.cpsize2(index_cloud_selected1(ii));
    p1.ctype(index_cloud_selected1(ii))   = 101;
    p1.ctype2(index_cloud_selected1(ii))  = -9999;    
    
    %% fix the cloud fracs
    p1.cfrac12(index_cloud_selected1(ii)) = 0;
    p1.cfrac2(index_cloud_selected1(ii))  = 0;
    p1.cngwat2(index_cloud_selected1(ii))  = 0;
    p1.cfrac(index_cloud_selected1(ii))   = 1;
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  for ii = 1 : length(index_cloud_selected3)
    p1.cfrac(index_cloud_selected3(ii))   = p1.cfrac2(index_cloud_selected3(ii));
    p1.cngwat(index_cloud_selected3(ii))  = p1.cngwat2(index_cloud_selected3(ii));
    p1.cprtop(index_cloud_selected3(ii))  = p1.cprtop2(index_cloud_selected3(ii));
    p1.cprbot(index_cloud_selected3(ii))  = p1.cprbot2(index_cloud_selected3(ii));
    p1.cpsize(index_cloud_selected3(ii))  = p1.cpsize2(index_cloud_selected3(ii));
    p1.ctype(index_cloud_selected3(ii))   = 101;
    p1.ctype2(index_cloud_selected3(ii))  = -9999;        
    
    %% fix the cloud fracs
    p1.cfrac12(index_cloud_selected3(ii)) = 0;
    p1.cfrac2(index_cloud_selected3(ii))  = 0;
    p1.cngwat2(index_cloud_selected3(ii))  = 0;
    p1.cfrac(index_cloud_selected3(ii))   = 1;
  end

  index_cloud_selected5 = find(p1.ctype  == 101 & p1.ctype2 == 101); %% so 1,2 are both water
  for ii = 1 : length(index_cloud_selected5)
    p1.clwc(:,index_cloud_selected5(ii)) = 0;
  
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(index_cloud_selected5(ii));
    cB = p1.cprbot(index_cloud_selected5(ii));
    xx = find(p1.plevs(:,index_cloud_selected5(ii)) >= cT & p1.plevs(:,index_cloud_selected5(ii)) <= cB);
    p1.clwc(xx,index_cloud_selected5(ii)) = sumwater(index_cloud_selected5(ii))/(cB-cT);
    p1.cfrac(index_cloud_selected5(ii))  = 1;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(index_cloud_selected5(ii));
    cB = p1.cprbot2(index_cloud_selected5(ii));
    xx = find(p1.plevs(:,index_cloud_selected5(ii)) >= cT & p1.plevs(:,index_cloud_selected5(ii)) <= cB);
    p1.clwc(xx,index_cloud_selected5(ii)) = sumwater(index_cloud_selected5(ii))/(cB-cT);
    p1.cfrac2(index_cloud_selected5(ii))  = 1;

    p1.cfrac12(index_cloud_selected5(ii))  = 1;
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waterORice < 0
  disp(' >>>>>>>>>>>>> SARTA : turning off water clouds <<<<<<<<<<<<<<')
  disp(' >>>> keeping only ice clouds ==> less profiles out than in <<<<')
  
  %% water is 101, ice is 201
  %% keep only ice slab clouds %%
  index_cloud_selected1 = find(p0.ctype  == 101 & p0.ctype2 == 201); %% so 1 = water, 2 = ice
  index_cloud_selected2 = find(p0.ctype2 == 101 & p0.ctype  == 201); %% so 1 = ice,   2 = water
  index_cloud_selected3 = find(p0.ctype2 == 201 & p0.ctype  <  101); %% so 2 = ice and 1 is not an water cloud (ie clear)
  index_cloud_selected4 = find(p0.ctype  == 201 & p0.ctype2 <  101); %% so 2 is not an water cloud (ie clear), 1 = ice
  index_cloud_selected5 = find(p0.ctype  == 201 & p0.ctype2 == 201); %% so 1,2 are both ice
  
  %% these are the ones we want yay yay yay
  index_cloud_selected = union(union(union(union(index_cloud_selected1,index_cloud_selected2),index_cloud_selected3),index_cloud_selected4),index_cloud_selected5);  

  %% these are the ones we do not want bah bah bah
  index_other = 1 : length(p0.stemp);
  index_other = setdiff(index_other,index_cloud_selected); %% either clear or both water or both ice

  junk = [length(index_cloud_selected1) length(index_cloud_selected2) length(index_cloud_selected3) length(index_cloud_selected4) length(index_cloud_selected5) length(index_other)];
  fprintf(1,'number of W/I I/W X/I I/X I/I X/X = %4i %4i %4i %4i %4i %4i \n',junk)

  [hx,p1x] = subset_rtp_allcloudfields(h,p0,[],[],index_cloud_selected);
  if length(index_other) > 0
    [hx,p1z] = subset_rtp_allcloudfields(h,p0,[],[],index_other);
    p1z.clwc = 0 * p1z.clwc;      %% set water to 0
    p1z.ciwc = 0 * p1z.ciwc;      %% set ice to 0 since it was not selected anyway!!! ???
  end
    
  p1 = p1x;
  p1.clwc = 0 * p1.clwc;        %% set water to 0
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumice(ii) = trapz(p1.plevs(1:nlevs,ii),p1.ciwc(1:nlevs,ii));
  end

  index_cloud_selected1 = find(p1.ctype  == 101 & p1.ctype2 == 201); %% so 1 = water, 2 = ice
  index_cloud_selected2 = find(p1.ctype2 == 101 & p1.ctype  == 201); %% so 1 = ice,   2 = water

  for ii = 1 : length(index_cloud_selected1)
    p1.cfrac(index_cloud_selected1(ii))   = 0;
    p1.cfrac12(index_cloud_selected1(ii)) = 0;
    p1.cngwat(index_cloud_selected1(ii))  = 0;
    p1.ctype(index_cloud_selected1(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(index_cloud_selected1(ii));
    cB = p1.cprbot2(index_cloud_selected1(ii));
    p1.ciwc(:,index_cloud_selected1(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected1(ii)) >= cT & p1.plevs(:,index_cloud_selected1(ii)) <= cB);
    p1.ciwc(xx,index_cloud_selected1(ii)) = sumice(index_cloud_selected1(ii))/(cB-cT);
  end
  
  for ii = 1 : length(index_cloud_selected2)
    p1.cfrac2(index_cloud_selected2(ii))  = 0;
    p1.cfrac12(index_cloud_selected2(ii)) = 0;
    p1.cngwat2(index_cloud_selected2(ii)) = 0;
    p1.ctype2(index_cloud_selected2(ii))  = -9999;
    
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(index_cloud_selected2(ii));
    cB = p1.cprbot(index_cloud_selected2(ii));
    p1.ciwc(:,index_cloud_selected2(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected2(ii)) >= cT & p1.plevs(:,index_cloud_selected2(ii)) <= cB);
    p1.ciwc(xx,index_cloud_selected2(ii)) = sumice(index_cloud_selected2(ii))/(cB-cT);
    p1.cfrac(index_cloud_selected2(ii))  = 1;
  end

  index_cloud_selected3 = find(p1.ctype2 == 201 & p1.ctype  < 101); %% so 2 = ice and 1 is not an water cloud (ie clear)
  index_cloud_selected4 = find(p1.ctype  == 201 & p1.ctype2 < 101); %% so 2 is not an water cloud (ie clear), 1 = ice

  for ii = 1 : length(index_cloud_selected3)
    p1.cfrac(index_cloud_selected3(ii))   = 0;
    p1.cfrac12(index_cloud_selected3(ii)) = 0;
    p1.cngwat(index_cloud_selected3(ii))  = 0;
    p1.ctype(index_cloud_selected3(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(index_cloud_selected3(ii));
    cB = p1.cprbot2(index_cloud_selected3(ii));
    p1.ciwc(:,index_cloud_selected3(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected3(ii)) >= cT & p1.plevs(:,index_cloud_selected3(ii)) <= cB);
    p1.ciwc(xx,index_cloud_selected3(ii)) = sumice(index_cloud_selected3(ii))/(cB-cT);
  end
  
  for ii = 1 : length(index_cloud_selected4)
    p1.cfrac2(index_cloud_selected4(ii))  = 0;
    p1.cfrac12(index_cloud_selected4(ii)) = 0;
    p1.cngwat2(index_cloud_selected4(ii)) = 0;
    p1.ctype2(index_cloud_selected4(ii))  = -9999;
    
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(index_cloud_selected4(ii));
    cB = p1.cprbot(index_cloud_selected4(ii));
    p1.ciwc(:,index_cloud_selected4(ii)) = 0;
    xx = find(p1.plevs(:,index_cloud_selected4(ii)) >= cT & p1.plevs(:,index_cloud_selected4(ii)) <= cB);
    p1.ciwc(xx,index_cloud_selected4(ii)) = sumice(index_cloud_selected4(ii))/(cB-cT);
    p1.cfrac(index_cloud_selected4(ii))  = 1;
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  for ii = 1 : length(index_cloud_selected1)
    p1.cfrac(index_cloud_selected1(ii))   = p1.cfrac2(index_cloud_selected1(ii));
    p1.cngwat(index_cloud_selected1(ii))  = p1.cngwat2(index_cloud_selected1(ii));
    p1.cprtop(index_cloud_selected1(ii))  = p1.cprtop2(index_cloud_selected1(ii));
    p1.cprbot(index_cloud_selected1(ii))  = p1.cprbot2(index_cloud_selected1(ii));
    p1.cpsize(index_cloud_selected1(ii))  = p1.cpsize2(index_cloud_selected1(ii));
    p1.ctype(index_cloud_selected1(ii))   = 201;
    p1.ctype2(index_cloud_selected1(ii))  = -9999;
    
    %% fix the cloud fracs
    p1.cfrac12(index_cloud_selected1(ii)) = 0;
    p1.cfrac2(index_cloud_selected1(ii))  = 0;
    p1.cngwat2(index_cloud_selected1(ii)) = 0;
    p1.cfrac(index_cloud_selected1(ii))   = 1;
  end

  for ii = 1 : length(index_cloud_selected3)
    p1.cfrac(index_cloud_selected3(ii))   = p1.cfrac2(index_cloud_selected3(ii));
    p1.cngwat(index_cloud_selected3(ii))  = p1.cngwat2(index_cloud_selected3(ii));
    p1.cprtop(index_cloud_selected3(ii))  = p1.cprtop2(index_cloud_selected3(ii));
    p1.cprbot(index_cloud_selected3(ii))  = p1.cprbot2(index_cloud_selected3(ii));
    p1.cpsize(index_cloud_selected3(ii))  = p1.cpsize2(index_cloud_selected3(ii));
    p1.ctype(index_cloud_selected3(ii))   = 201;
    p1.ctype2(index_cloud_selected3(ii))  = -9999;
    
    %% fix the cloud fracs
    p1.cfrac12(index_cloud_selected3(ii)) = 0;
    p1.cfrac2(index_cloud_selected3(ii))  = 0;
    p1.cngwat2(index_cloud_selected3(ii)) = 0;
    p1.cfrac(index_cloud_selected3(ii))   = 1;
  end

  index_cloud_selected5 = find(p1.ctype  == 201 & p1.ctype2 == 201); %% so 1,2 are both ice
  for ii = 1 : length(index_cloud_selected5)
    p1.ciwc(:,index_cloud_selected5(ii)) = 0;
  
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(index_cloud_selected5(ii))
    cB = p1.cprbot(index_cloud_selected5(ii))
    xx = find(p1.plevs(:,index_cloud_selected5(ii)) >= cT & p1.plevs(:,index_cloud_selected5(ii)) <= cB);
    p1.ciwc(xx,index_cloud_selected5(ii)) = sumice(index_cloud_selected5(ii))/(cB-cT);
    p1.cfrac(index_cloud_selected5(ii))  = 1;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(index_cloud_selected5(ii))
    cB = p1.cprbot2(index_cloud_selected5(ii))
    xx = find(p1.plevs(:,index_cloud_selected5(ii)) >= cT & p1.plevs(:,index_cloud_selected5(ii)) <= cB);
    p1.ciwc(xx,index_cloud_selected5(ii)) = sumice(index_cloud_selected5(ii))/(cB-cT);
    p1.cfrac2(index_cloud_selected5(ii))  = 1;

    p1.cfrac12(index_cloud_selected5(ii))  = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now merge
if length(index_other) > 0
  %% turn off any non-ice (waterORice = -1) or non-water (waterORice = +1) cloud
  p1z.cngwat  = 0 * ones(size(p1z.stemp));
  p1z.ctype   = -9999 * ones(size(p1z.stemp));
  p1z.cngwat2 = 0 * ones(size(p1z.stemp));
  p1z.ctype2  = -9999 * ones(size(p1z.stemp));
  
  [h1,p1] = merge_sort_rtp(hx,p1,hx,p1z,index_cloud_selected,index_other);
end

