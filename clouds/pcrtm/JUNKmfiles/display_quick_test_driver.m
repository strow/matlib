ii= 8;  ii= 35; 

ii = input('Enter ii : ');
[max(tmpjunk.ucol_num) tmpjunk.ucol_num(ii)]
%figure(1); semilogx(aa.p.ciwc(:,ii),aa.p.plevs(:,ii),'b+-',aa.p.clwc(:,ii),aa.p.plevs(:,ii),'ro-'); set(gca,'ydir','reverse')
%figure(2); pcolor(1:ncol0,aa.p.plevs(:,ii),squeeze(tmpjunk.unique_col_frac(ii,:,:))'); set(gca,'ydir','reverse'); colorbar
%figure(3); pcolor(1:ncol0,aa.p.plevs(:,ii),squeeze(tmpjunk.subcol_frac(ii,:,:))'); set(gca,'ydir','reverse'); colorbar

figure(1); semilogx(p.ciwc(:,ii),p.plevs(:,ii),'b+-',p.clwc(:,ii),p.plevs(:,ii),'ro-'); set(gca,'ydir','reverse')
  title('b = ice, r= water')
figure(2); plot(p.cc(:,ii),p.plevs(:,ii),'b+-');; set(gca,'ydir','reverse'); title('cloud frac')
figure(3); pcolor(1:ncol0,p.plevs(:,ii),squeeze(tmpjunk.unique_col_frac(ii,:,:))'); set(gca,'ydir','reverse'); colorbar
  title('unique col frac')
figure(4); pcolor(1:ncol0,p.plevs(:,ii),squeeze(tmpjunk.subcol_frac(ii,:,:))'); set(gca,'ydir','reverse'); colorbar
  title('subcol frac')
figure(5)
  pcolor(tmpjunk.ucol_num_same'); shading interp; colorbar; hold on; plot(tmpjunk.ucol_num,'k','linewidth',2); hold off
  tmpjunk.ucol_num_same(ii,:)