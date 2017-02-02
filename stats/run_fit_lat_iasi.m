fin = ['Data/Night/statlat'];
fout_pre = ['Data/fit_8year_night/fit_'];

%i = str2num(getenv('SLURM_ARRAY_TASK_ID'));

start_time = datenum(2007,09,24);
stop_time = datenum(2015,09,23);

parfor i=1:40
   fit_type = 'robs';
   fout = [fout_pre fit_type '_lat'];
   fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time);

   fit_type = 'rcal';
   fout = [fout_pre  'rcal_lat'];
   fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time);
end

% In a rush to see above, do bias last
parfor i=1:40
   fit_type = 'bias';
   fout = [fout_pre 'bias_lat'];
   fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time);
end
