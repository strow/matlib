%% Edit these two lines

fin = 'Data/Desc/statlat';
fout_pre = 'Data/Desc_fits/fit_';
%fout_pre = 'Data/Desc_fits/fit_5yr_a2c_';
% fin = 'Data/Desc/statlat';
% fout_pre = 'Data/Desc_fits/fit_';
% fin = 'Data/Desc_10K/statlat';
% fout_pre = 'Data/fit_13year_desc_10K/fit_';
%% Edit these two lines

start_time = datenum(2002,9,1,0,0,0);
stop_time = datenum(2018,8,31,0,0,0);
% start_time = datenum(2007,05,01);   % add 3 months and can include IASI!!
% stop_time = datenum(2017,05,01);
% % ---------------------------------------------------------------------------------- 
% fit_type = 'rcal';
% fout = [fout_pre 'rcal' '_lat'];
% for i=1:40
%    fit_robust_one_lat_stemp_subset(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
%---------------------------------------------------------------------------------- 
% % Mixed fit, 5 years a2c, 5 years cris
fit_type = 'bias';
fout = [fout_pre 'bias'  '_lat'];
for i=1:40
   fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'al1c');
   i
end
% % ---------------------------------------------------------------------------------- 
% % % For pure a2c fits, no CrIS involved
% fit_type = 'robs';
% fout = [fout_pre 'a2crobs' '_lat'];
% for i=22
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'a2cx');
%    i
% end
% % ---------------------------------------------------------------------------------- 
% fit_type = 'robs';
% fout = [fout_pre 'robs' '_lat'];
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

