%% Edit these two lines

fin = 'Data/Desc/statlat';
% fout_pre = 'Data/Asc_fits/fit_';
% fin = ['/r1/Data/Work/Airs/Random/Data/Asc/statlat'];
% fincal = ['/r1/Data/Work/Airs/Random/Data/Asc/nucal_statlat'];
fincal = [''];
fout_pre = 'Data/Desc_fits_cris_times/';
fout_presm = 'Data/Desc_fits_cris_times/Small/';
% fout_pre = ['Data/Desc_fits_iasi_times/fit_'];
% fout_presm = 'Data/Desc_fits_iasi_times/Small/fit_';

% fin = ['~/Work/Cris/Stability/Data/Desc/statlat'];
% fincal = [''];
% fout_pre = 'Data/Desc_fits_airs_times/';
% fout_presm = 'Data/Desc_fits_airs_times/Small/';


% fin = 'Data/Desc_10K/statlat';
% fout_pre = 'Data/fit_13year_desc_10K/fit_';
%% Edit these two lines

% % AIRS Good
% start_time = datenum(2002,9,1,0,0,0);
% stop_time = datenum(2018,08,31);

% Raw CrIS clear times
% start_time = datenum(2012,5,1,0,0,0);
% stop_time = datenum(2019,08,30,0,0,0);

% CrIS clear times to match AIRS (temp)
start_time = datenum(2012,5,1,0,0,0);
stop_time = datenum(2018,8,31,0,0,0);

% stop_time = datenum(2018,8,31,0,0,0);
%start_time = datenum(2007,09,24);


% start_time = datenum(2002,09,01);
% stop_time = datenum(2017,09,01);
% % ---------------------------------------------------------------------------------- 
fit_type = 'robs';
for i=1:40
   fout = [fout_pre 'robs_lat'];
   foutsm = [fout_presm 'robs' '_lat'];
%   fit_robust_one_lat_stemp_subset(fin,fout,i,fit_type,start_time,stop_time);
   fit_robust_one_lat(fin,fincal,fout,foutsm,i,fit_type,start_time,stop_time,'al1c','saveboth');
   i
end
% % ---------------------------------------------------------------------------------- 
% ---------------------------------------------------------------------------------- 
fit_type = 'rclr';
for i=1:40
   fout = [fout_pre 'rclr_lat'];
   foutsm = [fout_presm 'rclr' '_lat'];
%   fit_robust_one_lat_stemp_subset(fin,fout,i,fit_type,start_time,stop_time);
   fit_robust_one_lat(fin,fincal,fout,foutsm,i,fit_type,start_time,stop_time,'al1c','saveboth');
   i
end
% % ---------------------------------------------------------------------------------- 

% % ---------------------------------------------------------------------------------- 
fit_type = 'bias';
fout = [fout_pre 'bias' '_lat'];
foutsm = [fout_presm 'bias' '_lat'];
for i=1:40
   fit_robust_one_lat(fin,fincal,fout,foutsm,i,fit_type,start_time,stop_time,'al1c','saveboth');
   i
end
% ----------------------------------------------------------------------------------  


% fit_type = 'rclr';
% fout = [fout_pre 'rclr_resid' '_lat'];
% for i=1:40
% %   fit_robust_one_lat_stemp_subset(fin,fout,i,fit_type,start_time,stop_time);
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'al1c');
%    i
% end
% % % ---------------------------------------------------------------------------------- 
% fit_type = 'rcal';
% fout = [fout_pre 'rcal' '_lat'];
% for i=1:40
%    fit_robust_one_lat_stemp_subset(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% % ---------------------------------------------------------------------------------- 
% fit_type = 'bias';
% fout = [fout_pre 'bias' '_lat'];
% for i=1:40
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'al1c');
%    i
% end
% % ---------------------------------------------------------------------------------- 
% fit_type = 'rcld';
% fout = [fout_pre 'rcld' '_lat'];
% for i=1:40
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'al1c');
%    i
% end
% % ---------------------------------------------------------------------------------- 
% fit_type = 'rclr';
% fout = [fout_pre 'rclr' '_lat'];
% for i=1:40
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'al1c');
%    i
% end
% % % ---------------------------------------------------------------------------------- 
% fit_type = 'bias';
% fout = [fout_pre 'bias' '_lat'];
% for i=1:40
%    fit_robust_one_lat_subset(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% % ----------------------------------------------------------------------------------  
% fit_type = 'rcal';
% fout = [fout_pre 'rcal' '_lat'];
% for i=1:40
%    fit_robust_one_lat_subset(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% % ---------------------------------------------------------------------------------- 
% fit_type = 'rclr';
% fout = [fout_pre 'rclr' '_lat'];
% for i=1:40
%    fit_robust_one_lat_subset(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% ----------------------------------------------------------------------------------
% fit_type = 'bias';
% fout = [fout_pre 'cldforcing' '_lat'];
% for i=1:40
%    fit_robust_one_lat_cldforcing(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% 
% fit_type = 'rcal';
% fout = [fout_pre 'rcal' '_lat'];
% for i=1:40
%    fit_robust_one_lat_cldforcing(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% 
% 
% fit_type = 'rclr';
% fout = [fout_pre 'rclr' '_lat'];
% for i=1:40
%    fit_robust_one_lat_cldforcing(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end

