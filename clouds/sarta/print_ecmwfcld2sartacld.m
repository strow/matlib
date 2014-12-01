fprintf(1,' \n');
if (length(cTYPE) >= 1)
  disp('type    toplev   botlev  pTop pBot  cngwat*1000 kg/kg     cngwat g/m2');
  disp('-------------------------------------------------------------');
  for kk = 1 : length(cTYPE)
    if cTYPE(kk) == 'I'
      %fprintf(1,'  I      %3i       %3i         %8.6f         %8.6f \n',cT(kk),cB(kk),cngwat(kk)*1000,cc(kk))
      fprintf(1,'  I      %3i       %3i   %8.6f     %8.6f       %8.6f         %8.6f \n',cT(kk),cB(kk),plevs(cT(kk)),plevs(cB(kk)),sum(cOUT)*1000,cngwat(kk))
    else
      %fprintf(1,'  W      %3i       %3i         %8.6f         %8.6f \n',cT(kk),cB(kk),cngwat(kk)*1000,cc(kk))
      fprintf(1,'  W      %3i       %3i   %8.6f     %8.6f       %8.6f         %8.6f \n',cT(kk),cB(kk),plevs(cT(kk)),plevs(cB(kk)),sum(cOUT)*1000,cngwat(kk))
      if plevs(cT(kk)) < 420
        %error('lkjkj')
plot(profIN.ciwc(:,52),profIN.plevs(:,52),profIN.clwc(:,52),profIN.plevs(:,52),'r'); set(gca,'ydir','reverse'); grid; ax = axis; axis([0 ax(2) 0 1000])
        disp('ret to continue'); pause
      end
    end  %% if
  end    %% for
  disp('-------------------------------------------------------------');
else
  fprintf(1,' no cld found \n');
end     %% if length(ctype)
fprintf(1,' \n');



