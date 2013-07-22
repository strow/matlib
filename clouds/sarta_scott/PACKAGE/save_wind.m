iYear = num2str(xx(1));
iMonth = xx(2);
if iMonth < 10
  iMonth = ['0' num2str(iMonth)];
else
  iMonth = [    num2str(iMonth)];
  end
iDay = xx(3);
if iDay < 10
  iDay = ['0' num2str(iDay)];
else
  iDay = [    num2str(iDay)];
  end
iGran = xx(4);
if iGran < 10
  iGran = ['00' num2str(iGran)];
elseif iGran < 100
  iGran = ['0' num2str(iGran)];
else
  iGran = [    num2str(iGran)];
  end

saver = ['save winds_' iYear '_' iMonth '_' iDay '_' iGran ' XX YY GRAD ST U0 V0 U3 V3'];
eval(saver)

figure(1); caxis([270 310]); colorbar
figure(2); caxis([0 5]);     colorbar