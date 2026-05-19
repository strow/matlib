function [prof,hlayers,players] = main_compute_sarta_rads(h,ha,prof0,pa,pINPUT,run_sarta)

prof = prof0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  adds CO2 profile, if needed removes ice or water clds, runs sarta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if running interactively, and had old slabs but new slabs are computed,
%% can plot histograms
iCompareSlabs = +1;
iCompareSlabs = -1;
if iCompareSlabs > 0
  do_compare_plot_oldVSnew_slabs(prof,pINPUT,iCompareSlabs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add on co2
prof = prof_add_co2(h,prof,run_sarta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now if wanted can turn off ice or water clouds!!!!!!!!
%% this is kinda doing what "driver_sarta_cloud_rtp_onecldtest.m" is meant to do
if run_sarta.waterORice ~= 0
  disp('run_sarta.waterORice ~= 0 ---->>> going to remove ice or water cld')
  [prof,index_kept] = only_waterORice_cloud(h,prof,run_sarta.waterORice);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('inside driver_sarta_cloud_rtp.m --> main_compute_sarta_rads.m');
% which get_sarta_clear
% which get_sarta_cloud

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  %printarray([min(prof.rlon) max(prof.rlon) min(prof.rlat) max(prof.rlat)],'in main_compute_sarta_rads.m : min/max rlon  min.max rlat')
  [prof,hlayers,players] = get_sarta_clear(h,ha,prof,pa,run_sarta);
else
  disp('you did not ask for SARTA clear to be run; not changing prof.sarta_rclearcalc')  
end

if run_sarta.cloud > 0 
  disp('running SARTA cloud')
  %[prof.ctype  prof.cngwat  prof.cpsize  prof.cprtop  prof.cprbot]
  %[prof.ctype2 prof.cngwat2 prof.cpsize2 prof.cprtop2 prof.cprbot2]
  %figure(1); plot(prof.ciwc,prof.plevs);  ax = axis; line([ax(1) ax(2)],[prof.cprtop prof.cprtop],'color','k');   line([ax(1) ax(2)],[prof.cprbot prof.cprbot],'color','k');
  %figure(2); plot(prof.clwc,prof.plevs);  ax = axis; line([ax(1) ax(2)],[prof.cprtop2 prof.cprtop2],'color','k'); line([ax(1) ax(2)],[prof.cprbot2 prof.cprbot2],'color','k');
  [prof,hlayers,players] = get_sarta_cloud(h,ha,prof,pa,run_sarta);
else
  disp('you did not ask for SARTA cloudy to be run; not changing prof.rcalc')
end
