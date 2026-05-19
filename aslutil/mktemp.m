function tname = mktemp(TempDir,namePrefix);
%MKTEMP Create a unique temporary file name.
%   NAME = MKTEMP(dir,prefix), MKTEMP(dir_w_prefix), or MKTEMP()
%   MKTEMP uses the built-in unix system binary that generates 
%   unique temporary file names and verifies they don't
%   cause collisions.
%
%   Examples:
%      file1 = mktemp()
%      file2 = mktemp('/tmp/testme')
%      file3 = mktemp('/tmp','testme')
%      tmp_dir = mktemp('dir')  % an empty directory
%      ls /tmp
%      mktemp('clean')
%
%   Note: when specifying the directories /dev/shm this function will automatically
%     postfix the directory with a random number (per JOB / day span) on the cluster 
%     to avoid collisions.  The same will happen when /tmp is specified.
%
%   Warning: In MATLAB versions 2008a and later all files generated using 
%   this function will automatically be deleted when MatLab closes.

% Version 1.1, Written by Paul Schou - 26 Oct 2008
%  updated:  18 June 2011 - Paul Schou  added the ability to clear all files ending with name.*

persistent AlternateDir
global mktemp_handles

% cleanup routine for removing temporary files
if(nargin == 1 && strcmpi(TempDir,'clean'))
  if(verLessThan('matlab','7.6.0'))
    for i=1:length(mktemp_handles)
      disp(['  cleaning ' mktemp_handles{i}])
      unlink(mktemp_handles{i});
      for f = findfiles([mktemp_handles{i} '.*'])
        disp(['  cleaning extra ' f{1}])
        unlink(f{1});
      end
    end
    mktemp_handles = {};
  else
    mktemp_handles = [];
  end

  % Exit status if requested
  if(nargout == 1)
    tname=0;
  end
  return
elseif(nargin > 0 & strcmpi(TempDir,'dir'))
  if(nargin == 1)
    tname = mktemp();
  else
    tname = mktemp(namePrefix);
  end
  delete(tname);
  mkdir(tname);
  return;
end

if(nargout == 0)
  error('MKTEMP: no output variable specified');
end

%%% Main function %%%

  % handle the input arguments and assign variables
  if(nargin == 0)
    TempDir = getenv('JOB_SCRATCH_DIR');
    namePrefix = 'mktemp';
  elseif(nargin == 1)
    namePrefix = TempDir;
    TempDir = dirname(TempDir);
    if strcmp(TempDir,'.')
      TempDir = getenv('TMPDIR');          %% orig code
      TempDir = getenv('JOB_SCRATCH_DIR'); %% newer code      
      if(isdir('/scratch'))                %%% new Jan 2026
        [status, cmdout] = system('echo $SLURM_JOB_ID');
        if ~isempty(cmdout) && cmdout(end) == char(10)
          cmdout(end) = [];
        end
        TempDir = ['/scratch/' cmdout '/'];
      end      
    end
  elseif(~ isdir(TempDir) && isempty(AlternateDir))
    disp(['Warning: ' TempDir ' does not exist, using alternate temporary directory'])
    AlternateDir = 1; % supress future error messages
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%{
 Hi Sergio, typically as part of the job submission, scratch directories are
made automatically. See below:

Last login: Fri Jul 11 19:38:46 2025 from 10.2.49.213
[root@c24-52 ~]# cd /scratch/
[root@c24-52 scratch]# ls
234466 236800 236801 ebuild
[root@c24-52 scratch]#

[root@chip ~]# squeue | grep c24-52
234466 pi_strow interact sergio R 2-20:02:34 1 c24-52

Note that one of the directories has the same name as the job that's currently
running.

If I start a new job on the pi_strow a new directory will appear that
corresponds to the job that I'm running.

[sergio@chip ~]$ salloc --cluster=chip-cpu --account=pi_strow
--partition=pi_strow --qos=pi_strow --time=5 --mem=5G
salloc: Granted job allocation 236802
salloc: Waiting for resource configuration
salloc: Nodes c24-52 are ready for job

[sergio@c24-52 ~]$ cd /scratch
[sergio@c24-52 scratch]$ ls
234466 236801 236802 ebuild
[sergio@c24-52 scratch]$

A new directory was created with the same name as my job number. This is where
I can use the scratch storage. /tmp is a relatively small piece of storage in
comparison to the /scratch directory, and filling it can have unexpected
ramifications, so it's highly advised you not use it when possible. As part of
the epilogue of job, your scratch directories will automatically be cleaned up
when the job is cancelled or completed. In general though you can't create a
new directory in /scratch on your own. Hope this helps!

use      echo $SLURM_JOB_ID
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  % verify that the directory given is a valid location
  if(~ isdir(TempDir)) % if the directory does not exist, choose another
    if(isdir('/scratch'))    %%% new Jan 2026
      [status, cmdout] = system('echo $SLURM_JOB_ID');
      if ~isempty(cmdout) && cmdout(end) == char(10)
        cmdout(end) = [];
      end
      TempDir = ['/scratch/' cmdout '/'];            
    elseif(isdir('/tmp'))
      TempDir = '/tmp';
    elseif(isdir('/dev/shm'))
      TempDir = '/dev/shm';
    else
      error('MKTEMP:  No temporary folder found on system')
    end
  end

  if strcmp(TempDir(1:min(end,8)),'/dev/shm') & ~isempty(getenv('SHMDIR'))
    TempDir = getenv('SHMDIR');
  elseif strcmp(TempDir,'/tmp') & ~isempty(getenv('TMPDIR'))
    TempDir = getenv('TMPDIR');
  end
  
  % generate a random file name and verify it doesn't already exist
  validchars = [(0:25)+'a' (0:25)+'A' (0:9)+'0'];
  fid = fopen('/dev/urandom','r');
  while ( 1 )
    t=mod(fread(fid,8),62)+1;
    randval = validchars(t);
    [pathstr name ext] = fileparts(namePrefix);
    tname = [TempDir '/' name '_' randval ext];
    if ( ~ exist(tname,'file') )
      tfid = fopen(tname,'w+');
      if tfid
        fclose(tfid);
        break;
      else
        error(['MKTEMP: Could not create ' tname]);
      end
    end
  end
  fclose(fid);

  % add the file name to a cell array or the handle
  if(verLessThan('matlab','7.6.0'))
    mktemp_handles = {mktemp_handles , tname};
  else
    c = onCleanup(@()unlink(tname));
    mktemp_handles = [mktemp_handles c];
  end
