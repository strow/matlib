fin = ['Data/Desc/statlat'];
fout_pre = ['Data/Desc_fits/fit_'];

%i = str2num(getenv('SLURM_ARRAY_TASK_ID'));

start_time = datenum(2007,09,24);
stop_time = datenum(2018,09,31);

for i=2:40
   i
   fit_type = 'robs';
   fout = [fout_pre fit_type '_lat'];
   fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'iasi');

%    fit_type = 'rcld';
%    fout = [fout_pre  'rcld_lat'];
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'iasi');
end

% % In a rush to see above, do bias last
% for i=1:40
%    fit_type = 'bias';
%    fout = [fout_pre 'bias_lat'];
%    fit_robust_one_lat(fin,fout,i,fit_type,start_time,stop_time,'iasi');
% end
