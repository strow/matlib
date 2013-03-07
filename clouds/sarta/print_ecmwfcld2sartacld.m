fprintf(1,' \n');
if (length(cTYPE) >= 1)
  disp('type    toplev   botlev    cngwat*1000 kg/kg     cngwat g/m2');
  disp('-------------------------------------------------------------');
  for kk = 1 : length(cTYPE)
    if cTYPE(kk) == 'I'
      fprintf(1,'  I      %3i       %3i         %8.6f         %8.6f \n',cT(kk),cB(kk),cngwat(kk)*1000,cc(kk))
    else
      fprintf(1,'  W      %3i       %3i         %8.6f         %8.6f \n',cT(kk),cB(kk),cngwat(kk)*1000,cc(kk))
    end  %% if
  end    %% for
  disp('-------------------------------------------------------------');
else
  fprintf(1,' no cld found \n');
end     %% if length(ctype)
fprintf(1,' \n');
