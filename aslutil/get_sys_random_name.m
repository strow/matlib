function randname = get_sys_random_name();

% Schou: generate a random file name and verify it doesn't already exist
validchars = [(0:25)+'a' (0:25)+'A' (0:9)+'0'];
fid = fopen('/dev/urandom','r');
t=mod(fread(fid,8),62)+1;
randval = char(validchars(t));
fclose(fid);
tdir = getenv('JOB_SCRATCH_DIR');
randname = fullfile(tdir,randval);