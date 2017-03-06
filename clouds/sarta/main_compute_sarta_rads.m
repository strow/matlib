%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  adds CO2 profile, if needed removes ice or water clds, runs sarta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if running interactively, and had old slabs but new slabs are computed,
%% can plot histograms
do_compare_plot_oldVSnew_slabs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add on co2
prof_add_co2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now if wanted can turn off ice or water clouds!!!!!!!!
%% this is kinda doing what "driver_sarta_cloud_rtp_onecldtest.m" is meant to do
if run_sarta.waterORice ~= 0
  disp('run_sarta.waterORice ~= 0 ---->>> going to remove ice or water cld')
  [prof,index_kept] = only_waterORice_cloud(h,prof,run_sarta.waterORice);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_sarta.clear > 0 
  disp('running SARTA clear, saving into rclearcalc')
  tic
  get_sarta_clear;
  toc
  prof.sarta_rclearcalc = profRX2.rcalc;
else
  disp('you did not ask for SARTA clear to be run; not changing prof.sarta_rclearcalc')  
end

if run_sarta.cloud > 0 
  disp('running SARTA cloud')
  tic
  get_sarta_cloud;
  toc
  prof.rcalc = profRX2.rcalc;
else
  disp('you did not ask for SARTA cloudy to be run; not changing prof.rcalc')
end
