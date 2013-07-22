%% when you run readecmwf91_nearest_gasNcloud_slabprof in interactive mode

for ii = 1 : length(profX.cfrac)
  cprof = profX.cc(:,ii);
  clprof = profX.clwc(:,ii);
  ciprof = profX.ciwc(:,ii);
  cfrac1(ii) = sum(cprof.*clprof)/(sum(clprof)+1e-12);
  cfrac2(ii) = sum(cprof.*ciprof)/(sum(ciprof)+1e-12);

  wgtwater = profX.plevs(:,ii)/max(profX.plevs(:,ii));
  wgtwater = wgtwater.^4;

  wgtwater = zeros(length(profX.plevs(:,ii)),1);
    h1 = find(profX.plevs(:,ii) >= 600);
    wgtwater(h1) = 1;
    h2 = find(profX.plevs(:,ii) >= 400 & profX.plevs(:,ii) < 600);
    wgtwater(h2) = 0.66;
    h3 = find(profX.plevs(:,ii) >= 200 & profX.plevs(:,ii) < 400);
    wgtwater(h3) = 0.33;
  wgtice   = 1-wgtwater;

  cprof = profX.cc(:,ii);     %cprof = cprof.*cprof;
  clprof = profX.clwc(:,ii);  %clprof = clprof.*clprof;
  ciprof = profX.ciwc(:,ii);  %ciprof = ciprof.*ciprof;
  cfrac1(ii) = sum(cprof.*clprof.*wgtwater)/(sum(clprof.*wgtwater)+1e-12);
  cfrac2(ii) = sum(cprof.*ciprof.*wgtice)/(sum(ciprof.*wgtice)+1e-12);

  N = 10;
  cprof = profX.cc(:,ii);     cprof = cprof.^N;
  clprof = profX.clwc(:,ii);  %clprof = clprof.*clprof;
  ciprof = profX.ciwc(:,ii);  %ciprof = ciprof.*ciprof;
  cfrac1(ii) = (sum(cprof.*clprof)/(sum(clprof)+1e-12)).^(1/N);
  cfrac2(ii) = (sum(cprof.*ciprof)/(sum(ciprof)+1e-12)).^(1/N);


  h11 = subplot(131);
    plot([profX.cc(:,ii) wgtwater wgtice],profX.plevs(:,ii))
    set(gca,'ydir','reverse');
    axis([0 1 0 1200]);
    title(['cf' num2str(ii) ' = ' num2str(profX.cfrac(ii))]);
  h12 = subplot(132);
    plot(profX.clwc(:,ii),profX.plevs(:,ii))
    set(gca,'ydir','reverse');
    axis([0 max(profX.clwc(:,ii))*1.001+1e-8 0 1200]);
    title(['W']);
  h13 = subplot(133);
    plot(profX.ciwc(:,ii),profX.plevs(:,ii))
    set(gca,'ydir','reverse');
    axis([0 max(profX.ciwc(:,ii))*1.001+1e-8 0 1200]);
    title(['I']);
  [ii max(profX.cc(:,ii)) profX.cfrac(ii) cfrac1(ii) cfrac2(ii)]
  pause;
  end
