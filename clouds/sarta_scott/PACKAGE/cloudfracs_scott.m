cfrac = profX.cfrac;
cc    = profX.cc;
ciwc  = profX.ciwc;
clwc  = profX.clwc;
[ciout,cwout,fraci,fracw,fraciw] = fudge_cloud(cfrac,cc,ciwc,clwc);

for ii = 1 : length(profX.cfrac)
  h11 = subplot(131);
    plot(profX.cc(:,ii),profX.plevs(:,ii))
    set(gca,'ydir','reverse');
    axis([0 1 0 1200]);
    title(['cf' num2str(ii) ' = ' num2str(profX.cfrac(ii))]);
  h12 = subplot(132);
    plot(profX.clwc(:,ii),profX.plevs(:,ii),cwout(:,ii),profX.plevs(:,ii))
    set(gca,'ydir','reverse');
    axis([0 max(profX.clwc(:,ii))*1.001+1e-8 0 1200]);
    title(['W']);
  h13 = subplot(133);
    plot(profX.ciwc(:,ii),profX.plevs(:,ii),ciout(:,ii),profX.plevs(:,ii))
    set(gca,'ydir','reverse');
    axis([0 max(profX.ciwc(:,ii))*1.001+1e-8 0 1200]);
    title(['I']);
  [max(profX.cc(:,ii)) profX.cfrac(ii) fracw(ii) fraci(ii)]
  pause;
  end


