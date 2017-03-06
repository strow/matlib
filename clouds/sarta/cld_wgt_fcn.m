function [wgtW,wgtpeakWindex,wgtpeakW_tempr,wgtpeakW,wgtI,wgtpeakIindex,wgtpeakI_tempr,wgtpeakI] = cld_wgt_fcn(p1);

wgtI = zeros(size(p1.ptemp));
wgtW = zeros(size(p1.ptemp));

%% do the weighting functions
ii = 1;
wgtI(ii,:) = 0;
for ii = 2 : p1.nlevs(1)
  wgtI(ii,:) = wgtI(ii-1,:) + p1.sarta_lvlODice(ii,:);
end
wgtI = exp(-wgtI) .* (1-exp(-p1.sarta_lvlODice));
for ii = 1 : length(p1.stemp)
  dodo = wgtI(:,ii);
  if sum(dodo) > 0
    bop = find(dodo == max(dodo),1);
    wgtpeakIindex(ii) = bop;
    wgtpeakI_tempr(ii) = p1.ptemp(bop,ii);
    wgtpeakI(ii) = p1.plevs(bop,ii);
  else
    wgtpeakIindex(ii) = -1;
    wgtpeakI_tempr(ii) = -9999;
    wgtpeakI(ii) = -9999;
  end
end

ii = 1;
wgtI(ii,:) = 0;
for ii = 2 : p1.nlevs(1)
  wgtW(ii,:) = wgtW(ii-1,:) + p1.sarta_lvlODwater(ii,:);
end
wgtW = exp(-wgtW) .* (1-exp(-p1.sarta_lvlODwater));
for ii = 1 : length(p1.stemp)
  dodo = wgtW(:,ii);
  if sum(dodo) > 0
    bop = find(dodo == max(dodo),1);
    wgtpeakWindex(ii) = bop;
    wgtpeakW_tempr(ii) = p1.ptemp(bop,ii);
    wgtpeakW(ii) = p1.plevs(bop,ii);
  else
    wgtpeakWindex(ii) = -1;
    wgtpeakW_tempr(ii) = -9999;
    wgtpeakW(ii) = -9999;
  end
end
