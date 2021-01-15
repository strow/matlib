fin = ['Data/Desc/statlat'];
fout_pre = ['Data/Desc_fits/fit_'];
fout_presm = 'Data/Desc_fits/Small/fit_';

%i = str2num(getenv('SLURM_ARRAY_TASK_ID'));

start_time = datenum(2007,09,24);
stop_time = datenum(2018,09,31);

% for i=1:40
%    i
%    fit_type = 'robs';
%    fout = [fout_pre fit_type '_lat'];
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'iasi');

%    fit_type = 'rcld';
%    fout = [fout_pre  'rcld_lat'];
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'iasi');
% end

% In a rush to see above, do bias last
for ifov = 1:4
for i=4:35
   fit_type = 'robs';
   fout = [fout_pre 'robs_lat_' int2str(ifov) '_'];
   foutsm = [fout_presm 'robs' '_lat_' int2str(ifov) '_'];
   fit_robust_one_lat(fin,fout,foutsm,i,fit_type,start_time,stop_time,'iasd','saveboth',ifov);
end
end


