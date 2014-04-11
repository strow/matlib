function [head, hattr, prof, pattr] = rtpklayers_keepfields(head, hattr, prof, pattr, klayers_exec)
%function [head, hattr, prof, pattr] = rtpklayers_keepfields(head, hattr, prof, pattr, klayers_exec)
%
%  A simple function to run klayers on a file and return the result
%  Exacutable is in pattr field "klayers_exec".
%
%  This routine keeps all provided PROF fields by running klayers only on the
%  profile relevant routines and replacing it in the original PROF.
%
%  klayers_exec (optional) name of the executable - overrides pattr field.
%
% P.S./B.I. - 2014/04/12


% Are we already layer profile?
if head.ptype > 0
  disp('  rtpklayers: Already in levels, doing nothing')
  return
end

% Separate all the profile relevant fields:
[h0 p0] = keep_klfields_l(head,prof);


% Run Klayers in the smaller RTP structure
if(nargin()<5)
  klayers_exec = get_attr(hattr,'klayers_exec');
end

if ~isempty(klayers_exec)

  tmp1 = mktemp();
  tmp2 = mktemp();
  rtpwrite(tmp1, h0, hattr, p0, pattr);

  disp(['    ' klayers_exec ' fin=' tmp1 ' fout=' tmp2 ' > /dev/null']);
  system([klayers_exec ' fin=' tmp1 ' fout=' tmp2 ' > /dev/null']);
  delete(tmp1)

  [h0, hattr, p0, pattr] = rtpread(tmp2);
  delete(tmp2)
else
  if(nargin()<5)
    error('Cannot run klayers as klayers_exec variable is missing from hattr.')
  else
    error('Cannot run klayers as klayers_exec input argument is empty.')
  end
end


% Now merge only the relevant fields
[head prof ] = merge_klfields_l(head, prof, h0, p0);

% Done!

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h1 p1] = keep_klfields_l(h,p)
 
  % list of relevant klayers fields - non-gases
  fn={'plat','plon','ptime','stemp','salti','spres','nlevs','plevs','ptemp','palts','gtotal','gxover','txover','co2ppm'};
  for iff=1:numel(fn)
    if(isfield(p,fn{iff}))
      p1.(fn{iff})=p.(fn{iff});
    end
  end

  % Do gases
  for ig=1:h.ngas
    gn=h.glist(ig);
    p1.(['gas_' num2str(gn)]) = p.(['gas_' num2str(gn)]);
  end

  % Fix header
  h1=h;
  h1.pfields=1;
  h1.nchan=0;
  if(isfield(h1,'ichan'))
    h1=rmfield(h1,'ichan');
  end
  if(isfield(h1,'vchan'))
    h1=rmfield(h1,'vchan');
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hbase pbase] = merge_klfields_l(hbase, pbase, hkl, pkl)

  % Will add the Layer relevant fields from hkl and pkl into
  % hbase and pbase

  fn={'plat','plon','ptime','stemp','salti','spres','nlevs','plevs','ptemp','palts','gtotal','gxover','txover','co2ppm'};
  for iff=1:numel(fn)
    if(isfield(pkl,fn{iff}))
      pbase.(fn{iff})=pkl.(fn{iff});
    end
  end

  % Do gases
  for ig=1:hkl.ngas
    gn=hkl.glist(ig);
    pbase.(['gas_' num2str(gn)]) = pkl.(['gas_' num2str(gn)]);
  end


  % Fix header
  hbase.ptype = hkl.ptype;
  hbase.ngas = hkl.ngas;
  hbase.glist = hkl.glist;
  hbase.gunit = hkl.gunit;

end 
  




