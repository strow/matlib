thedir = dir('*.m');
thedirX = dir('/home/sergio/MATLABCODE/matlib/clouds/sarta/*.m');

if length(thedir) ~= length(thedirX)
  error('oops different number of files ....')
end

for ii = 1 : length(thedir)
  fname1 = thedir(ii).name;
  fname2 = ['/home/sergio/MATLABCODE/matlib/clouds/sarta/' fname1];
  ee1 = exist(fname1);
  ee2 = exist(fname2);
  if ee2 == 0
    fprintf(1,'>>>>>>>>> %s does not exist in the matlib dir\n',fname1);
  else
    fprintf(1,'file %3i out of %3i : diffing %s \n',ii,length(thedir),fname1);
    differer = ['!diff ' fname1 ' ' fname2];
    eval(differer);
    disp('ret to continue');
    pause
  end
end
