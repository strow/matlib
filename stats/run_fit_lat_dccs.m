%% Edit these two lines

% To do nucal, only explicitely set in fout_pre and in call (fit_robust*)
fin = 'Data/Desc/statlat';
%fout_pre = 'Data/Desc_fits/fit_nucal_';
fout_pre = 'Data/Desc_fits_iasi_times/fit_';
% fin = 'Data/Desc_merra/statlat';
% fout_pre = 'Data/Desc_merra_fits/fit_';

% fin = 'Data/Desc_10K/statlat';
% fout_pre = 'Data/fit_13year_desc_10K/fit_';
%% Edit these two lines

% start_time = datenum(2012,4,24,0,0,0);
% stop_time = datenum(2016,4,23,0,0,0);
%start_time = datenum(2002,09,01);
start_time = datenum(2007,07,24);
stop_time = datenum(2016,07,23);

fit_type = 'robs';
fout = [fout_pre 'robs' '_lat'];
for i=14:27
   fit_robust_one_lat_dccs_airs(fin,fout,i,fit_type,start_time,stop_time);
%   fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time);
   i
end
% 
% fit_type = 'rcal';
% fout = [fout_pre 'rcal' '_lat'];
% for i=1:40
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% 
% % fit_type = 'cf';
% % fout = [fout_pre 'cf' '_lat'];
% % for i=1:40
% % fit_robust_one_lat_nucal_cf_cal(fin,fout,i,fit_type,start_time,stop_time);
% % %   fit_robust_one_lat_nucal
% %    i
% % end
% % 
% fit_type = 'rclr';
% fout = [fout_pre 'rclr' '_lat'];
% for i=1:40
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% 
% fit_type = 'bias';
% fout = [fout_pre 'bias' '_lat'];
% for i=1:40
%    fit_robust_one_lat_nucal(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end

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
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time);
%    i
% end
% 
