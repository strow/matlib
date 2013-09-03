function [head, hattr, prof, pattr] = rtpklayers(head, hattr, prof, pattr, klayers_exec)
%function [head, hattr, prof, pattr] = rtpklayers(head, hattr, prof, pattr, klayers_exec)
%
%  A simple function to run klayers on a file and return the result
%  Exacutable is in pattr field "klayers_exec".
%
% klayers_exec (optional) name of the executable - overrides pattr field.

%if isfield(head,'gunit') && head.gunit(1) == 1
if head.ptype > 0
  disp('  rtpklayers: Already in levels, doing nothing')
  return
end

if isfield(prof,'rcalc'); rcalc = prof.rcalc; prof = rmfield(prof,'rcalc'); end
if isfield(prof,'robs1'); robs1 = prof.robs1; prof = rmfield(prof,'robs1'); end
if isfield(prof,'calflag'); calflag = prof.calflag; prof = rmfield(prof,'calflag'); end
h = head;
head.pfields = 1;
head.nchan = 1;
head.vchan = 1;
head.ichan = 1;

if(nargin()<5)
  klayers_exec = get_attr(hattr,'klayers_exec');
end

if ~isempty(klayers_exec)
  tmp1 = mktemp();
  tmp2 = mktemp();
  rtpwrite(tmp1, head, hattr, prof, pattr);
  clear head hattr prof pattr
  disp(['    ' klayers_exec ' fin=' tmp1 ' fout=' tmp2 ' > /dev/null']);
  system([klayers_exec ' fin=' tmp1 ' fout=' tmp2 ' > /dev/null']);
  delete(tmp1)
  [head, hattr, prof, pattr] = rtpread(tmp2);
  delete(tmp2)
else
  if(nargin()<5)
    error('Cannot run klayers as klayers_exec variable is missing from hattr.')
  else
    error('Cannot run klayers as klayers_exec input argument is empty.')
  end
end

head.pfields = h.pfields;
head.nchan = h.nchan;
head.ichan = h.ichan;
head.vchan = h.vchan;
if exist('robs1','var'); prof.robs1 = robs1; end
if exist('rcalc','var'); prof.rcalc = rcalc; end
if exist('calflag','var'); prof.calflag = calflag; end
