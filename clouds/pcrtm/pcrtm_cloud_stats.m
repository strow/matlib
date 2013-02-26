%disp('doing the PCRTM cloud stats ....')

clear totalODwater meanDMEwater maxCTOPwater da_w
clear totalODice meanDMEice maxCTOPice da_i

for icol = 2:ucol_num(ibox)
  clear da_w da_i

  da_w = find(cldphase(icol,:) == 1); 
  if length(da_w) > 1
    totalODwater(icol) = sum(cldopt(icol,da_w));
    meanDMEwater(icol) = mean(cldde(icol,da_w));
    maxCTOPwater(icol) = min(cldpres(icol,da_w));
  elseif length(da_w) == 1
    totalODwater(icol) = (cldopt(icol,da_w));
    meanDMEwater(icol) = (cldde(icol,da_w));
    maxCTOPwater(icol) = (cldpres(icol,da_w));
  else
    totalODwater(icol) = NaN;
    meanDMEwater(icol) = NaN;
    maxCTOPwater(icol) = NaN;
  end

  da_i = find(cldphase(icol,:) == 2); 
  if length(da_i) > 1
    totalODice(icol) = sum(cldopt(icol,da_i));
    meanDMEice(icol) = mean(cldde(icol,da_i));
    maxCTOPice(icol) = min(cldpres(icol,da_i));
  elseif length(da_i) == 1
    totalODice(icol) = (cldopt(icol,da_i));
    meanDMEice(icol) = (cldde(icol,da_i));
    maxCTOPice(icol) = (cldpres(icol,da_i));
  else
    totalODice(icol) = NaN;
    meanDMEice(icol) = NaN;
    maxCTOPice(icol) = NaN;
  end
end

%tmpjunk
%whos total* meanDME* maxCTOP*
%  [ibox ucol_num(ibox)]

%% assume no ice or water in any of the Ncol0 profiles
tmpjunk.totalODice(ibox) = NaN;
tmpjunk.meanDMEice(ibox) = NaN;
tmpjunk.maxCTOPice(ibox) = NaN;
tmpjunk.totalODwater(ibox) = NaN;
tmpjunk.meanDMEwater(ibox) = NaN;
tmpjunk.maxCTOPwater(ibox) = NaN;

%% now check for ice or water
if ucol_num(ibox) > 1
  ixy = 2:ucol_num(ibox);

  isi = find(isfinite(meanDMEice(ixy)));
  if length(isi) > 0
    summ = sum(ucol_num_same(ibox,ixy(isi)));
    tmpjunk.totalODice(ibox) = nansum(totalODice(ixy(isi)).*ucol_num_same(ibox,ixy(isi)))/summ;
    tmpjunk.meanDMEice(ibox) = nansum(meanDMEice(ixy(isi)).*ucol_num_same(ibox,ixy(isi)))/summ;
    tmpjunk.maxCTOPice(ibox) = nansum(maxCTOPice(ixy(isi)).*ucol_num_same(ibox,ixy(isi)))/summ;
  end

  isw = find(isfinite(meanDMEwater(ixy)));
  if length(isw) > 0
    summ = sum(ucol_num_same(ibox,ixy(isw)));
    tmpjunk.totalODwater(ibox) = nansum(totalODwater(ixy(isw)).*ucol_num_same(ibox,ixy(isw)))/summ;
    tmpjunk.meanDMEwater(ibox) = nansum(meanDMEwater(ixy(isw)).*ucol_num_same(ibox,ixy(isw)))/summ;
    tmpjunk.maxCTOPwater(ibox) = nansum(maxCTOPwater(ixy(isw)).*ucol_num_same(ibox,ixy(isw)))/summ;
  end
end
