function [head hattr prof pattr] = rtpadd_merra_interp(head,hattr,prof,pattr)
% function [head hattr prof pattr] = rtpadd_merra(head,hattr,prof,pattr)


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fetch MERRA Data 
  % Fields 	description	RTP field name
  % t  - air temperature  	ptemp
  % qv - specific humidity	gas_1
  % o3 - ozone mixing ratio	gas_3
  % ps - surface pressure		spres
  % ts - surface temperature	stemp
  % u2m - eastward wind at 2 m above the displacement height
  % v2m - northward wind at 2 m above the displacement height
  %
  % Data is given on 3hr files, except 'ts' which is hourly
  %


  % Assumptions:
  % 1. All pressure grids are the same for all 3D variables
  % 2. Invalid bottom of atmosphere values happens on all 3D variables


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Basic checks
  if(head.ptype~=0)
    disp('You have an RTP file that is a level file. All previous layer information will be removed');
    while head.ngas>0
      if(isfield(prof,['gas_' num2str(head.glist(1))]))
	prof=rmfield(prof,['gas_' num2str(head.glist(1))]);
	head.ngas  = head.ngas-1;
	head.glist = head.glist(2:end,1);
	head.gunit = head.gunit(2:end,1);
	head.ptype = 0;
      else
	warning(['Non existing field gas_' num2str(head.glist(1)) ' indicated by headers. Fixing']);
	head.ngas  = head.ngas-1;
	head.glist = head.glist(2:end,1);
	head.gunit = head.gunit(2:end,1);
	head.ptype = 0;
      end
    end
  end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get the data times:
  mtimes = AirsDate(prof.rtime); % [mtimes]=days

  nprof=numel(mtimes); 
  threehours = round((mtimes-mtimes(1))*8); % [threehours]=3-hour long units
  frac_hour=mtimes-floor(mtimes(1));   % fraction of the day  
  frac_interval=frac_hour*24;   % Fraction of the day in hours?  
  daily_time_1=floor(frac_hour*24/3)*3;  % hour of the day for preceding analysis time
  daily_time_2=floor(frac_hour*24/3)*3 +3; % hour of the day for following analysis time
  u_daily_time_1=unique(daily_time_1); 
  u_daily_time_2=unique(daily_time_2); 
  time_1=floor(frac_hour*24/3)*3/24+floor(mtimes(1)); % preceding analysis time for each profile 
  time_2=time_1+3/24;
%  time_2=ceil(frac_hour*24/3)*3+floor(mtimes(nprof));  % following analysis time for each profile
  u_time_1=unique(time_1);  % unique list of the first analsys time
  u_time_2=unique(time_2);  % unique list of the 2nd analysis time
  n3hours=numel(u_time_1);  % Number of different time intervals
  frac_anal=(frac_interval-daily_time_1)/3; % Fractional time of each profile between the 2 analysis times 


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Loop over 
  
  for i3hours=1:n3hours
  i3hours

    % Select Fovs and subset them (using first analysis time)
    ifov=find(time_1==u_time_1(i3hours));   % Profiles with preceding analysis time of u_time_1
    nfovs = numel(ifov);
    for i=1:nfovs 
      frac_anal_new(i)=frac_anal(ifov(i)); 
    end 
    tprof = ProfSubset2(prof,ifov);

    % Get required analysis times before and after profile times 
%    reqtime_1 = floor(mtimes(1)) + u_time_1(i3hours)/24; % [day]=[day]+[Analysis time]
   
    reqtime_1=u_time_1(i3hours); 
    reqtime_2=u_time_2(i3hours); 
    'reqtime_1=',reqtime_1
    'reqtime_2=',reqtime_2 
%    reqtime=mtimes(1)+u_time_1(i3hours)/8;
%    reqtime_2 = floor(mtimes(1)) + u_time_2(i3hours)/24; % same  


    %%%%%%%%%%%%%%%%%%%% 
    % Set profile variables

    % ptime
    tprof.ptime_1 = ones(1,nfovs).*AirsDate(reqtime_1,-1); % [sec]
    tprof.ptime_2 = ones(1,nfovs).*AirsDate(reqtime_2,-1); 
    tprof.plat = tprof.rlat;
    tprof.plon = tprof.rlon;

    ptemp=[]; pgas_1=[]; pgas_3=[]; pps=[]; pts=[];
    pgas_1_t1=[]; pgas_1_t2=[]; pgas_3_t1=[]; pgas_3_t2=[]; 
    w2m=[]; 


    %%%%%%%%%%%%%%%%%%%% 
    % Interpolate 3D variables for each layer
    % ATTENTION: Fill value is not consistent. Some times it is 1e15 (as
    % advertised, but some times it is -9.99e8. Go figure!
 
    % t  - air temperature  	ptemp
    [dat_t_1 plevs_1 lats lons]= getdata_merra(reqtime_1, 't');
    [dat_t_2 plevs_2 lats lons]= getdata_merra(reqtime_2, 't'); 

    nlevs=numel(plevs_1);

    dat_t_1(dat_t_1>1e14 | dat_t_1<-1)=NaN;
    dat_t_2(dat_t_2>1e14 | dat_t_2<-1)=NaN; 
    for ilev=1:nlevs
      ptemp_1(ilev,:) = interp2(lats, lons, dat_t_1(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
      ptemp_2(ilev,:) = interp2(lats, lons, dat_t_2(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest'); 
      for itemp=1:numel(ifov);
          iprof=ifov(itemp); % ptemp access the subset array, but ptemp_1 and ptemp_2 need to start from 1. 
          ptemp(ilev,itemp) = ptemp_1(ilev,itemp) + (ptemp_2(ilev,itemp) - ptemp_1(ilev,itemp))*frac_anal(iprof); 
      end 
    end
    clear ptemp_1
    clear ptemp_2 
    
     
    
    % qv - specific humidity	gas_1
    'interp','i3hours=',i3hours 
    [dat_q_t1 plevs_1 lats lons]= getdata_merra(reqtime_1, 'qv');
    [dat_q_t2 plevs_2 lats lons]= getdata_merra(reqtime_2, 'qv'); 
    dat_q_t1(dat_q_t1>1e14 | dat_t_1<-1)=NaN;
    dat_q_t2(dat_q_t2>1e14 | dat_t_1<-1)=NaN;
    for ilev=1:nlevs
      pgas_1_t1(ilev,:) = interp2(lats, lons, dat_q_t1(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
      pgas_1_t2(ilev,:) = interp2(lats, lons, dat_q_t2(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest'); 
      for itemp=1:numel(ifov);
           iprof=ifov(itemp);
           pgas_1(ilev,itemp) = pgas_1_t1(ilev,itemp) + (pgas_1_t2(ilev,itemp) - pgas_1_t1(ilev,itemp))*frac_anal(iprof);
      end
    end
    clear pgas_1_t1;
    clear pgas_1_t2; 
    
   
    % o3 - ozone mixing ratio	gas_3
    [dat_o3_t1 plevs_1 lats lons]= getdata_merra(reqtime_1, 'o3');
    [dat_o3_t2 plevs_2 lats lons]= getdata_merra(reqtime_2, 'o3'); 
    dat_o3_t1(dat_o3_t1>1e14 | dat_o3_t1 < -1 | dat_t_1<-1)=NaN;
    dat_o3_t2(dat_o3_t2>1e14 | dat_o3_t1 < -1 | dat_t_2<-1)=NaN; 
    for ilev=1:nlevs
      pgas_3_t1(ilev,:) = interp2(lats, lons, dat_o3_t1(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
      pgas_3_t2(ilev,:) = interp2(lats, lons, dat_o3_t2(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest'); 
      for itemp=1:numel(ifov);
          iprof=ifov(itemp); 
          pgas_3(ilev,itemp) = pgas_3_t1(ilev,itemp) + (pgas_3_t2(ilev,itemp) - pgas_3_t1(ilev,itemp))*frac_anal(iprof);
      end
    end

    plevs_1=plevs_1(end:-1:1);
    plevs_2=plevs_2(end:-1:1); 
    plevs=plevs_1(end:-1:1);    % For now, use plevs from first analysis time, should interpolate eventually  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Compute Valid Levels
   
    % For each profile, find the last valid level. 
    % 
    tprof.plevs = plevs_1*ones(1, nfovs);
    tprof.nlevs = nlevs.*zeros(1, nfovs);
    for itemp=1:numel(ifov)
      iprof=itemp;
      tprof.nlevs(1,iprof) = nlevs-sum(isnan(ptemp(:,iprof)));
    end 

    tprof.ptemp = ptemp;
    tprof.gas_1 = pgas_1;
    tprof.gas_3 = pgas_3;
    


    % ps - surface pressure	spres
    [dat_ps_1 tmpx lats lons]= getdata_merra(reqtime_1, 'ps'); % It is in Pa, convert to mbar -> /100
    [dat_ps_2 tmpx lats lons]= getdata_merra(reqtime_2, 'ps'); 
    dat_ps_1(dat_ps_1>1e14)=NaN;
    dat_ps_2(dat_ps_2>1e14)=NaN; 
    pps_1 = interp2(lats,lons,dat_ps_1/100, tprof.rlat, tprof.rlon, 'nearest');
    pps_2 = interp2(lats,lons,dat_ps_2/100, tprof.rlat, tprof.rlon, 'nearest'); 
    for itemp=1:numel(ifov);
          iprof=ifov(itemp);  
          pps(itemp) = pps_1(itemp) + (pps_2(itemp) - pps_1(itemp))*frac_anal(iprof);
    end

    
    % ts - surface temperature	stemp
    [dat_ts_1 tmpx lats lons merra_str]= getdata_merra(reqtime_1, 'ts');
    [dat_ts_2 tmpx lats lons merra_str]= getdata_merra(reqtime_2, 'ts'); 
    dat_ts_1(dat_ts_1>1e14)=NaN;
    dat_ts_2(dat_ts_2>1e14)=NaN; 
    pts_1 = interp2(lats,lons,dat_ts_1, tprof.rlat, tprof.rlon, 'nearest');
    pts_2 = interp2(lats,lons,dat_ts_2, tprof.rlat, tprof.rlon, 'nearest'); 
    for itemp=1:numel(ifov);
           iprof=ifov(itemp); 
           pts(itemp) = pts_1(itemp) + (pts_2(itemp) - pts_1(itemp))*frac_anal(iprof);
    end
   
    tprof.spres = pps;
    tprof.stemp = pts;  

    
    % wind speed at 2m
    [dat_u2m_1 tmpx lats lons]= getdata_merra(reqtime_1, 'u2m');
    [dat_u2m_2 tmpx lats lons]= getdata_merra(reqtime_2, 'u2m'); 
    [dat_v2m_1 tmpx lats lons]= getdata_merra(reqtime_1, 'v2m');
    [dat_v2m_2 tmpx lats lons]= getdata_merra(reqtime_2, 'v2m'); 
    dat_w2m_1 = sqrt(dat_u2m_1.^2 + dat_v2m_1.^2);
    dat_w2m_2 = sqrt(dat_u2m_2.^2 + dat_v2m_2.^2); 
    w2m_1 = interp2(lats,lons,dat_w2m_1,tprof.rlat, tprof.rlon,'nearest'); 
    w2m_2 = interp2(lats,lons,dat_w2m_2,tprof.rlat, tprof.rlon,'nearest'); 
    for itemp=1:numel(ifov);
        iprof=ifov(itemp); 
        w2m(itemp) = w2m_1(itemp) + (w2m_2(itemp) - w2m_1(itemp))*frac_anal(iprof);
    end 
    clear w2m_1;
    clear w2m_2; 

    tprof.wspeed = w2m;


    tprof_arr(i3hours) = tprof;

  end
  tprof = Prof_join_arr(tprof_arr);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fix Header and Attributes


  % Set that the profile type is level.
  head.ptype=0;

  % Set the profile bit
  [ia ib ic]=pfields2bits(head.pfields);
  head.pfields = bits2pfields(1,ib,ic);

  % Add the new two gases - 1 and 3
  for ig=[1 3]
    if(~isfield(head,'glist'))
      head.glist=[];
      head.gunit=[];
      head.ngas=0;
    end
    ik=find(head.glist==ig); 
    if(numel(ik)==0)
      head.glist(end+1,1)=ig; % gad id
      head.gunit(end+1,1)=21; % gas unit (g/g) or something like that
      head.ngas=head.ngas+1;
    else
      head.gunit(ik)=21;
    end 
  end

  % Set pressures
  head.pmin = min(tprof.plevs(:));
  head.pmax = max(tprof.spres(:));

  hattr = set_attr(hattr,'profile',[merra_str ' Nearest']);

  prof=tprof;

end
