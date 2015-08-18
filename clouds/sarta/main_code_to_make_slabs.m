%% first do some sanity checks

oo = find(p.cc < 0);   p.cc(oo) = 0.0;
oo = find(p.cc > 1);   p.cc(oo) = 1.0;
oo = find(p.ciwc < 0); p.ciwc(oo) = 0.0;
oo = find(p.clwc < 0); p.clwc(oo) = 0.0;

if sum(size(p.ptemp) - size(p.cc)) ~= 0 | sum(size(p.ptemp) - size(p.ciwc)) ~= 0 | sum(size(p.ptemp) - size(p.clwc)) ~= 0
  pjunk.cc    = p.cc;
  pjunk.ciwc  = p.ciwc;
  pjunk.clwc  = p.clwc;
  pjunk.ptemp = p.ptemp;
  pjunk
  error('sizes of ptemp cc ciwc clwc seem different')
end
if sum(size(p.ptemp) - size(p.gas_1)) ~= 0 | sum(size(p.ptemp) - size(p.gas_3)) ~= 0 | sum(size(p.ptemp) - size(p.plevs)) ~= 0
  pjunk.gas_1 = p.gas_1;
  pjunk.gas_3 = p.gas_3;
  pjunk.plevs = p.plevs;
  pjunk.ptemp = p.ptemp;
  pjunk  
  error('sizes of ptemp plevs gas_1 gas_3 seem different')
end
if isfield(p,'gas_2')
  if sum(size(p.gas_1) - size(p.gas_2)) ~= 0
    pjunk.gas_1 = p.gas_1;
    pjunk.gas_2 = p.gas_2;    
    pjunk  
    error('sizes of gas_1 gas_2 seem different')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlev     = ceil(mean(p.nlevs));
nlev_std = (std(double(p.nlevs)));

if h.ptype ~= 0
  error('need levels input!')
end

%if nlev_std > 1e-3
%  error('oops : code assumes ERA (37 levs) or ECMWF (91 levs) or other constant numlevs model')
%end

load airsheights.dat
load airslevels.dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

head  = h;

tic
if nlev_std <= 1e-3
  [prof,profX] = ecmwfcld2sartacld(p,nlev,run_sarta.cumsum,airslevels,airsheights);   %% figure the two slab cloud 
                    %% profile info here, using profX
                    %% this then puts the info into "prof" by calling put_into_prof w/in routine
else
  disp('oops : code assumes ERA (37 levs) or ECMWF (91 levs) or other constant numlevs model, need to use varying levels (MERRA??)')
  [prof,profX] = ecmwfcld2sartacld_varNlev(p,nlev,run_sarta.cumsum,airslevels,airsheights);   %% figure the two slab cloud 
                    %% profile info here, using profX
                    %% this then puts the info into "prof" by calling put_into_prof w/in routine
end

prof = put_into_V201cld_fields(prof);    %% puts cloud info from above into rtpv201 fields 
  prof.ctype  = double(prof.ctype);
  prof.ctype2 = double(prof.ctype2);

if iDebugMain > 0
  figure(3); plot(prof.ciwc,prof.plevs,'b',prof.clwc,prof.plevs,'r'); ax=axis; axis([0 ax(2) 0 1000]); set(gca,'ydir','reverse')
  figure(1); plot_hists_cprtop_cprtop2
  disp('new keyboard posn1')
  keyboard
end

%% sets fracs and particle effective sizes eg cfrac2
prof = set_fracs_deffs(head,prof,profX,cmin,cngwat_max,run_sarta.cfrac,run_sarta.randomCpsize);

if run_sarta.cumsum > 0 & run_sarta.cumsum <= 1
  %% set cloud top according to cumulative sum fraction of ciwc or clwc
  profXYZ = prof;
  prof = reset_cprtop(prof);
elseif run_sarta.cumsum > 1
  %% set cloud top according to where cumulative cloud OD = run_sarta.cumsum/100, 
  %%      or if >= 9999, set at peak of cloud wgt fcn
  profXYZ = prof;
  prof = reset_cprtop_cloudOD(prof,run_sarta.cumsum/100,airslevels,airsheights);  
end

if iDebugMain > 0
  figure(2); plot_hists_cprtop_cprtop2
  figure(4); plot(prof.sarta_wgtI,prof.plevs,'b',prof.sarta_wgtW,prof.plevs,'r'); ax=axis; axis([0 ax(2) 0 1000]); 
           set(gca,'ydir','reverse')
  figure(5); 
    plot(prof.sarta_index_wgtpeakI(oo11),prof.cprtop(oo11),'bo',prof.sarta_index_wgtpeakI(oo12),prof.cprtop2(oo12),'cx',...
         prof.sarta_index_wgtpeakI(oo21),prof.cprtop(oo21),'ro',prof.sarta_index_wgtpeakW(oo22),prof.cprtop2(oo22),'mx')
             set(gca,'ydir','reverse')
    xlabel('layer of peak wgtfcn'); ylabel('cprtop')
  disp('new keyboard posn2')
  keyboard
end

disp('---> checking cprtop vs cprbot vs spres')
iNotOK = +1;
iFix = 0;
%% before used to give up at iFix == 10
while iFix < 12 & iNotOK > 0
  iFix = iFix + 1;
  [prof,iNotOK] = check_for_errors(prof,run_sarta.cfrac,iFix);  %% see possible pitfalls in clouds
  fprintf(1,' did n=%2i try at checking clouds \n',iFix)
end
if iFix >= 12 & iNotOK > 0
  %disp('oops, could not fix cprtop vs cprbot vs spres'); %keyboard
  error('oops, could not fix cprtop vs cprbot vs spres')
end

clear profX

if run_sarta.cfrac >= 0 & run_sarta.cfrac <= 1
  oo = find(prof.cfrac  > 0 & prof.cngwat > 0);  prof.cfrac(oo)  = run_sarta.cfrac;
  oo = find(prof.cfrac2 > 0 & prof.cngwat2 > 0); prof.cfrac2(oo) = run_sarta.cfrac;
  oo = find(prof.cfrac  > 0 & prof.cfrac2 > 0);  prof.cfrac12(oo)  = run_sarta.cfrac;
end
