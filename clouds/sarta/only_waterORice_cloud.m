function [p1,select_clouds,select_other] = only_waterORice_cloud(h,p0,waterORice)

%% p0 is input profile
%%      "waterORice" +1 find pixels with (one water slab & one ice slab), removes the ice slab cloud, rejigs water profile
%%                      find pixels with (one water slab)                                           , rejigs water profile
%%                   -1 find pixels with (one water slab & one ice slab), removes the water slab cloud, rejigs ice profile
%%                      find pixels with (one ice slab)                                             , rejigs ice profile
%%                    0 do   nothing
%% p1 is output profile,
%%   select_clouds are those whose clouds were selected as the pixel had (one water and one ice slab)  only
%%   select_other  are those that are clear or both ice or both water
%% 
%% >>>>>>>>>>>
%% so basically if you want to test PCRTM vs SARTA SLAB cloud you should
%{
[h,ha,p,pa] = colocate_era_airs(yy,mm,dd,gg); or make_rtp_yada(yy,mm,dd,gg);
run_sarta.ncol0 = -1;      %%% instead of 50, use -1 for testing ONE cloud
run_sarta.waterORice = +1; %%% turns off ice, keeps water clouds
run_sarta.waterORice = -1; %%% turns off water, keeps ice clouds

%% sub-select according to water or ice; re-jigs relevant cloud profile
[p1,select_clouds,select_other] = only_waterORice_cloud(h,p0,run_sarta.waterORice);

%% call sarta and pcrtm on the SLAB profiles !!!
p2 = driver_pcrtm_cloud_rtp(h,ha,p1,pa,run_sarta);                                        
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waterORice == 0
  disp(' >>>>>>>>>>>>> in only_waterORice_cloud.m, doing NOTHING (keeping ice/water clds) <<<<<<<<<<<<<<<<')
  p1 = p0;
  select_clouds = [];
  select_other  = [];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(p0.ctype == 101) == 0 & length(p0.ctype == 201) == 0 & length(p0.ctype2 == 101) == 0 & length(p0.ctype2 == 201) == 0
  error('only_waterORice_cloud : assumes profiles already mapped to slabs, now can select the I or W clouds and change profiles')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load airslevels.dat
load airsheights.dat
aLH.airslevels  = airslevels;
aLH.airsheights = airsheights;

if waterORice > 0
  disp(' >>>>>>> SARTA : keep water clds, turning off ice clds <<<<<<<<<')
  
  %% water is 101, ice is 201
  %% keep only water slab clouds %%
  select_clouds1 = find(p0.ctype  == 201 & p0.ctype2 == 101); %% so 2 = water, 1 = ice
  select_clouds2 = find(p0.ctype2 == 201 & p0.ctype  == 101); %% so 2 = ice,   1 = water
  select_clouds3 = find(p0.ctype2 == 101 & p0.ctype   < 101); %% so 2 = water and 1 is not an ice cloud (ie clear)
  select_clouds4 = find(p0.ctype  == 101 & p0.ctype2  < 101); %% so 2 is not an ice cloud (ie clear), 1 = water
  select_clouds5 = find(p0.ctype  == 101 & p0.ctype2 == 101); %% so 1,2 are both water

  %% these are the ones we want yay yay yay
  select_clouds = union(union(union(union(select_clouds1,select_clouds2),select_clouds3),select_clouds4),select_clouds5);

  %% these are the ones we do not want bah bah bah
  select_other = 1 : length(p0.stemp);
  select_other = setdiff(select_other,select_clouds); %% either clear or both water or both ice

  junk = [length(select_clouds1) length(select_clouds2) ...
          length(select_clouds3) length(select_clouds4) length(select_clouds5) length(select_other)];
  fprintf(1,'number of I/W W/I X/W W/X W/W  -- X/X = %4i %4i %4i %4i %4i -- %4i \n',junk)
  if sum(junk) ~= length(p0.rlat)
    fprintf(1,'sum(partition) : %5i total num profs : %5i \n',sum(junk),length(p0.rlat))
    error('sum of cloud partitions <> number of input profiles')
  end
  
  [hx,p1x] = subset_rtp_allcloudfields(h,p0,[],[],select_clouds);
  if length(select_other) > 0
    [hx,p1z] = subset_rtp_allcloudfields(h,p0,[],[],select_other);
    p1z.ciwc = 0 * p1z.ciwc;      %% set ice to 0
    p1z.clwc = 0 * p1z.clwc;      %% set water to 0, since it was not selected anyway!!! ???
  end
  
  p1 = p1x;
  p1.ciwc = 0 * p1.ciwc;        %% set ice to 0
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumwater(ii) = trapz(p1.plevs(1:nlevs,ii),p1.clwc(1:nlevs,ii));
  end

  select_clouds1 = find(p1.ctype  == 201 & p1.ctype2 == 101); %% so 2 = water, 1 = ice
  select_clouds2 = find(p1.ctype2 == 201 & p1.ctype  == 101); %% so 2 = ice,   1 = water
  select_clouds3 = find(p1.ctype2 == 101 & p1.ctype  < 101); %% so 2 = water and 1 is not an ice cloud (ie clear)
  select_clouds4 = find(p1.ctype  == 101 & p1.ctype2 < 101); %% so 2 is not an ice cloud (ie clear), 1 = water
  select_clouds5 = find(p1.ctype  == 101 & p1.ctype2 == 101); %% so 1,2 are both water

  indx = select_clouds1; 
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat(indx(ii))  = 0;
    p1.cprtop(indx(ii))  = -9999;    
    p1.cprbot(indx(ii))  = -9999;        
    p1.ctype(indx(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere    
    cT = p1.cprtop2(indx(ii));
    cB = p1.cprbot2(indx(ii));
    p1.cfrac2(indx(ii))   = 1;
    %% old
    p1.clwc(:,indx(ii)) = 0;    
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.clwc(xx,indx(ii)) = sumwater(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.clwc(:,indx(ii)) = 0;
    p1.clwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat2(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);
  end

  indx = select_clouds2;
  for ii = 1 : length(indx)
    p1.cfrac2(indx(ii))  = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat2(indx(ii)) = 0;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(indx(ii));
    cB = p1.cprbot(indx(ii));
    p1.cfrac(indx(ii))   = 1;        
    %% old    
    p1.clwc(:,indx(ii)) = 0;
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.clwc(xx,indx(ii)) = sumwater(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.clwc(:,indx(ii)) = 0;
    p1.clwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    
  end

  indx = select_clouds3;
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat(indx(ii))  = 0;
    p1.cprtop(indx(ii))  = -9999;    
    p1.cprbot(indx(ii))  = -9999;            
    p1.ctype(indx(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere    
    cT = p1.cprtop2(indx(ii));
    cB = p1.cprbot2(indx(ii));
    p1.cfrac2(indx(ii))  = 1;    
    %% old    
    p1.clwc(:,indx(ii)) = 0;
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.clwc(xx,indx(ii)) = sumwater(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.clwc(:,indx(ii)) = 0;
    p1.clwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat2(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    
  end

  indx = select_clouds4;
  for ii = 1 : length(indx)
    p1.cfrac2(indx(ii))  = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat2(indx(ii)) = 0;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(indx(ii));
    cB = p1.cprbot(indx(ii));
    p1.cfrac(indx(ii))  = 1;
    %% old    
    p1.clwc(:,indx(ii)) = 0;
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.clwc(xx,indx(ii)) = sumwater(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.clwc(:,indx(ii)) = 0;
    p1.clwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  indx = select_clouds1;
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = p1.cfrac2(indx(ii));
    p1.cngwat(indx(ii))  = p1.cngwat2(indx(ii));
    p1.cprtop(indx(ii))  = p1.cprtop2(indx(ii));
    p1.cprbot(indx(ii))  = p1.cprbot2(indx(ii));
    p1.cpsize(indx(ii))  = p1.cpsize2(indx(ii));
    p1.ctype(indx(ii))   = 101;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;    
    
    %% fix the cloud fracs
    p1.cfrac12(indx(ii)) = 0;
    p1.cfrac2(indx(ii))  = 0;
    p1.cngwat2(indx(ii))  = 0;
    p1.cfrac(indx(ii))   = 1;
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  indx = select_clouds3;  
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = p1.cfrac2(indx(ii));
    p1.cngwat(indx(ii))  = p1.cngwat2(indx(ii));
    p1.cprtop(indx(ii))  = p1.cprtop2(indx(ii));
    p1.cprbot(indx(ii))  = p1.cprbot2(indx(ii));
    p1.cpsize(indx(ii))  = p1.cpsize2(indx(ii));
    p1.ctype(indx(ii))   = 101;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;        
    
    %% fix the cloud fracs
    p1.cfrac12(indx(ii)) = 0;
    p1.cfrac2(indx(ii))  = 0;
    p1.cngwat2(indx(ii))  = 0;
    p1.cfrac(indx(ii))   = 1;
  end

  indx = select_clouds5;
  for ii = 1 : length(indx)
    p1.clwc(:,indx(ii)) = 0;
  
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(indx(ii));
    cB = p1.cprbot(indx(ii));
    p1.cfrac(indx(ii))  = 1;
    %% old    
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.clwc(xx,indx(ii)) = sumwater(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.clwc(:,indx(ii)) = 0;
    p1.clwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(indx(ii));
    cB = p1.cprbot2(indx(ii));
    p1.cfrac2(indx(ii))  = 1;    
    %% old    
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.clwc(xx,indx(ii)) = sumwater(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.clwc(:,indx(ii)) = 0;
    p1.clwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat2(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    

    p1.cfrac12(indx(ii))  = 1;
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waterORice < 0
  disp(' >>>>>>> SARTA : keep ice clds, turning off water clds <<<<<<<<<')
  
  %% water is 101, ice is 201
  %% keep only ice slab clouds %%
  select_clouds1 = find(p0.ctype  == 101 & p0.ctype2 == 201); %% so 1 = water, 2 = ice
  select_clouds2 = find(p0.ctype2 == 101 & p0.ctype  == 201); %% so 1 = ice,   2 = water
  select_clouds3 = find(p0.ctype2 == 201 & p0.ctype  <  101); %% so 2 = ice and 1 is not an water cloud (ie clear)
  select_clouds4 = find(p0.ctype  == 201 & p0.ctype2 <  101); %% so 2 is not an water cloud (ie clear), 1 = ice
  select_clouds5 = find(p0.ctype  == 201 & p0.ctype2 == 201); %% so 1,2 are both ice

  %% these are the ones we want yay yay yay
  select_clouds = union(union(union(union(select_clouds1,select_clouds2),select_clouds3),select_clouds4),select_clouds5);  

  %% these are the ones we do not want bah bah bah
  select_other = 1 : length(p0.stemp);
  select_other = setdiff(select_other,select_clouds); %% either clear or both water or both ice

  junk = [length(select_clouds1) length(select_clouds2) ...
          length(select_clouds3) length(select_clouds4) length(select_clouds5) length(select_other)];
  fprintf(1,'number of W/I I/W X/I I/X I/I -- X/X = %4i %4i %4i %4i %4i -- %4i \n',junk)
  if sum(junk) ~= length(p0.rlat)
    fprintf(1,'sum(partition) : %5i total num profs : %5i \n',sum(junk),length(p0.rlat))
    error('sum of cloud partitions <> number of input profiles')
  end
  
  [hx,p1x] = subset_rtp_allcloudfields(h,p0,[],[],select_clouds);
  if length(select_other) > 0
    [hx,p1z] = subset_rtp_allcloudfields(h,p0,[],[],select_other);
    p1z.clwc = 0 * p1z.clwc;      %% set water to 0
    p1z.ciwc = 0 * p1z.ciwc;      %% set ice to 0 since it was not selected anyway!!! ???
  end
    
  p1 = p1x;
  p1.clwc = 0 * p1.clwc;        %% set water to 0
  for ii = 1 : length(p1.stemp)
    nlevs = p1.nlevs(ii);
    sumice(ii) = trapz(p1.plevs(1:nlevs,ii),p1.ciwc(1:nlevs,ii));
  end

  select_clouds1 = find(p1.ctype  == 101 & p1.ctype2 == 201); %% so 1 = water, 2 = ice
  select_clouds2 = find(p1.ctype2 == 101 & p1.ctype  == 201); %% so 1 = ice,   2 = water
  select_clouds3 = find(p1.ctype2 == 201 & p1.ctype  < 101); %% so 2 = ice and 1 is not an water cloud (ie clear)
  select_clouds4 = find(p1.ctype  == 201 & p1.ctype2 < 101); %% so 2 is not an water cloud (ie clear), 1 = ice
  select_clouds5 = find(p1.ctype  == 201 & p1.ctype2 == 201); %% so 1,2 are both ice

  indx = select_clouds1;
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat(indx(ii))  = 0;
    p1.cprtop(indx(ii))  = -9999;    
    p1.cprbot(indx(ii))  = -9999;            
    p1.ctype(indx(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(indx(ii));
    cB = p1.cprbot2(indx(ii));
    p1.cfrac2(indx(ii))  = 1;
    %% old    
    p1.ciwc(:,indx(ii)) = 0;
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.ciwc(xx,indx(ii)) = sumice(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.ciwc(:,indx(ii)) = 0;
    p1.ciwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat2(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    
  end

  indx = select_clouds2;
  for ii = 1 : length(indx)
    p1.cfrac2(indx(ii))  = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat2(indx(ii)) = 0;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;
    
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(indx(ii));
    cB = p1.cprbot(indx(ii));
    p1.cfrac(indx(ii))  = 1;    
    %% old    
    p1.ciwc(:,indx(ii)) = 0;
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.ciwc(xx,indx(ii)) = sumice(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.ciwc(:,indx(ii)) = 0;
    p1.ciwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    
  end

  indx = select_clouds3;
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat(indx(ii))  = 0;
    p1.cprtop(indx(ii))  = -9999;    
    p1.cprbot(indx(ii))  = -9999;            
    p1.ctype(indx(ii))   = -9999;

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(indx(ii));
    cB = p1.cprbot2(indx(ii));
    p1.cfrac2(indx(ii))  = 1;    
    %% old    
    p1.ciwc(:,indx(ii)) = 0;
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.ciwc(xx,indx(ii)) = sumice(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.ciwc(:,indx(ii)) = 0;
    p1.ciwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat2(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    
  end

  indx = select_clouds4;
  for ii = 1 : length(indx)
    p1.cfrac2(indx(ii))  = 0;
    p1.cfrac12(indx(ii)) = 0;
    p1.cngwat2(indx(ii)) = 0;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;
    
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(indx(ii));
    cB = p1.cprbot(indx(ii));
    p1.cfrac(indx(ii))  = 1;    
    %% old    
    p1.ciwc(:,indx(ii)) = 0;
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.ciwc(xx,indx(ii)) = sumice(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.ciwc(:,indx(ii)) = 0;
    p1.ciwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);    
  end

  %% now make sure all info is in cngwat, cprtop etc and not in cngwat2 cprtop2 etc
  indx = select_clouds1;
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = p1.cfrac2(indx(ii));
    p1.cngwat(indx(ii))  = p1.cngwat2(indx(ii));
    p1.cprtop(indx(ii))  = p1.cprtop2(indx(ii));
    p1.cprbot(indx(ii))  = p1.cprbot2(indx(ii));
    p1.cpsize(indx(ii))  = p1.cpsize2(indx(ii));
    p1.ctype(indx(ii))   = 201;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;
    
    %% fix the cloud fracs
    p1.cfrac12(indx(ii)) = 0;
    p1.cfrac2(indx(ii))  = 0;
    p1.cngwat2(indx(ii)) = 0;
    p1.cfrac(indx(ii))   = 1;
  end

  indx = select_clouds3;
  for ii = 1 : length(indx)
    p1.cfrac(indx(ii))   = p1.cfrac2(indx(ii));
    p1.cngwat(indx(ii))  = p1.cngwat2(indx(ii));
    p1.cprtop(indx(ii))  = p1.cprtop2(indx(ii));
    p1.cprbot(indx(ii))  = p1.cprbot2(indx(ii));
    p1.cpsize(indx(ii))  = p1.cpsize2(indx(ii));
    p1.ctype(indx(ii))   = 201;
    p1.cprtop2(indx(ii)) = -9999;    
    p1.cprbot2(indx(ii)) = -9999;            
    p1.ctype2(indx(ii))  = -9999;
    
    %% fix the cloud fracs
    p1.cfrac12(indx(ii)) = 0;
    p1.cfrac2(indx(ii))  = 0;
    p1.cngwat2(indx(ii)) = 0;
    p1.cfrac(indx(ii))   = 1;
  end

  indx = select_clouds5;
  for ii = 1 : length(indx)
    p1.ciwc(:,indx(ii)) = 0;
  
    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop(indx(ii));
    cB = p1.cprbot(indx(ii));
    p1.cfrac(indx(ii))  = 1;    
    %% old    
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.ciwc(xx,indx(ii)) = sumice(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.ciwc(:,indx(ii)) = 0;
    p1.ciwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);

    %% try to change the input n-level profile into a slab profile, 0 elsewhere
    cT = p1.cprtop2(indx(ii));
    cB = p1.cprbot2(indx(ii));
    p1.cfrac2(indx(ii))  = 1;    
    %% old    
    xx = find(p1.plevs(:,indx(ii)) >= cT & p1.plevs(:,indx(ii)) <= cB);
    if length(xx) == 0
      xx = find(p1.plevs(:,indx(ii)) >= cB,1);
    end    
    p1.ciwc(xx,indx(ii)) = sumice(indx(ii))/(cB-cT);
    %% new
    pT = cT;
    pB = cB;
    p1.ciwc(:,indx(ii)) = 0;
    p1.ciwc(xx,indx(ii)) = convert_gm2_to_gg(pT,pB,p1.cngwat2(indx(ii)),...
                                              p1.plevs(:,indx(ii)),p1.ptemp(:,indx(ii)),p1.nlevs(indx(ii)),aLH);
    
    p1.cfrac12(indx(ii))  = 1;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now merge

if length(select_other) > 0
  %% turn off any non-ice (waterORice = -1) or non-water (waterORice = +1) cloud

  p1z.cfrac   = -9999 * ones(size(p1z.stemp));
  p1z.cngwat  = -9999 * ones(size(p1z.stemp));
  p1z.ctype   = -9999 * ones(size(p1z.stemp));
  p1z.cprtop  = -9999 * ones(size(p1z.stemp));
  p1z.cprbot  = -9999 * ones(size(p1z.stemp));

  p1z.cfrac2  = -9999 * ones(size(p1z.stemp));
  p1z.cngwat2 = -9999 * ones(size(p1z.stemp));
  p1z.ctype2  = -9999 * ones(size(p1z.stemp));
  p1z.cprtop2 = -9999 * ones(size(p1z.stemp));
  p1z.cprbot2 = -9999 * ones(size(p1z.stemp));

  p1z.cfrac12 = -9999 * ones(size(p1z.stemp));

  p1z.ciwc    = 0 * p1z.ciwc;
  p1z.clwc    = 0 * p1z.clwc;
  p1z.cc      = 0 * p1z.cc;
  
  [h1,p1] = merge_sort_rtp(hx,p1,hx,p1z,select_clouds,select_other);
end

