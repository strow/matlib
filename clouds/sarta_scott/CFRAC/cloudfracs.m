%% when you run readecmwf91_nearest_gasNcloud_slabprof in interactive mode

for ii = 1 : length(profX.cfrac)
  cfrac1(ii) = ...
    sum(profX.cc(:,ii).*profX.clwc(:,ii))/(sum(profX.clwc(:,ii))+1e-12);
  cfrac2(ii) = ...
    sum(profX.cc(:,ii).*profX.ciwc(:,ii))/(sum(profX.ciwc(:,ii))+1e-12);

  h11 = subplot(131);
    plot(profX.cc(:,ii),profX.plevs(:,ii))
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
  [max(profX.cc(:,ii)) profX.cfrac(ii) cfrac1(ii) cfrac2(ii)]
  pause;
  end


