%% so eg nboxes = 100 = number of FOVS to process (this is iChunk in driver_pcrtm_cloud_rtp.m)
%% while ibox = where in this loop are you????

%% also see get_subcolumn_frac_v2.m

%% >>>>>>>>> look at subr get_subcolumn_frac_v2.m <<<<<<<<<
% whos subcol_frac unique_col_frac ucol_num
jett = jet; jett(1,:) = 1;

fprintf(1,'looking at MRO clouds : index %4i of %4i has %3i unique sub_columns \n',ibox,nboxes,ucol_num(ibox))

figure(1); clf;
  lalalala = cc(:,ibox); lalalala = 1:length(lalalala);
  subplot(121); plot(ncol*cc(:,ibox),lalalala); title('NWP cc');
    ax1 = axis; axis([ax1(1) ax1(2) 1 length(lalalala)])
    set(gca,'ydir','reverse')
  subplot(122); plot(ICT(:,ibox),lalalala,'b',WCT(:,ibox),lalalala,'r'); hold off;
    title('NWP CIWC,CLWC')
    ax2 = axis; axis([ax2(1) ax2(2) 1 length(lalalala)])
    set(gca,'ydir','reverse')

%% do max, random, max random
[xunique_col_frac,xucol_num,xucol_num_same,xsubcol_frac] = get_subcolumn_frac_v2_debug(1, nlev, ncol, cc(:,ibox)', 1);
  figure(2); clf; imagesc(squeeze(xsubcol_frac(1,:,:))'); colormap(jett); colorbar; title('Maximum')
  axis([ax1(1) ax1(2) 1 length(lalalala)])
[xunique_col_frac,xucol_num,xucol_num_same,xsubcol_frac] = get_subcolumn_frac_v2_debug(1, nlev, ncol, cc(:,ibox)', 2);
  figure(3); clf; imagesc(squeeze(xsubcol_frac(1,:,:))'); colormap(jett); colorbar; title('Random')
[xunique_col_frac,xucol_num,xucol_num_same,xsubcol_frac] = get_subcolumn_frac_v2_debug(1, nlev, ncol, cc(:,ibox)', 3);
  figure(4); clf; imagesc(squeeze(xsubcol_frac(1,:,:))'); colormap(jett); colorbar; title('Maximum Random')

%whos xunique_col_frac xucol_num xucol_num_same xsubcol_frac
%sum(xucol_num_same)-ncol
%lalalala = cc(:,ibox); lalalala = 1:length(lalalala);
%  plot(squeeze(xsubcol_frac(1,1,:)),lalalala,cc(:,ibox),lalalala,'r'); axis([0 1 1 max(lalalala)])

figure(1); clf;
  lalalala = cc(:,ibox); lalalala = 1:length(lalalala);
  subplot(121); plot(cc(:,ibox),lalalala); title('NWP cc');
    ax = axis; axis([ax(1) ax(2) 1 length(lalalala)])
    set(gca,'ydir','reverse')
  subplot(122); plot(ICT(:,ibox),lalalala,'b',WCT(:,ibox),lalalala,'r'); hold off;
    title('NWP CIWC,CLWC')
    ax = axis; axis([ax(1) ax(2) 1 length(lalalala)])
    set(gca,'ydir','reverse')

if ucol_num(ibox) > 1
  figure(2); clf
    pcolor(squeeze(unique_col_frac(ibox,:,:))'); shading flat; colormap(jett); colorbar; title('cld frac')
    ylabel('vertical')
    set(gca,'ydir','reverse')
    
  figure(3); clf
    subplot(121); pcolor(cld_qw'); shading flat; colormap(jett); colorbar; title('cld qw')
    set(gca,'ydir','reverse')    
      ylabel('vertical')
    subplot(122); pcolor(cld_qi'); shading flat; colormap(jett); colorbar; title('cld qi')
    set(gca,'ydir','reverse')
    
  figure(4); clf
    subplot(121); pcolor(cldphase'); shading flat; colormap(jett); colorbar; title('cld phase')
    set(gca,'ydir','reverse')    
      ylabel('vertical')
      %% 1 = ice, 2 = water
    subplot(122); pcolor(cld_id'); shading flat; colormap(jett); colorbar; title('cld ID')
    set(gca,'ydir','reverse')
    
  figure(5); clf
    subplot(121); pcolor(cldopt'); shading flat; colormap(jett); colorbar; title('cld OD')
    set(gca,'ydir','reverse')    
      ylabel('vertical')
    subplot(122); pcolor(cldde'); shading flat; colormap(jett); colorbar; title('cld eff diam')
    set(gca,'ydir','reverse')
    
  figure(6); clf; pcolor(cldpres'); shading flat; colormap(jett); colorbar; title('cld pres')
    set(gca,'ydir','reverse')
    
  [junk,lalalala] = size(cldpres);

  figure(7); clf; pcolor(ones(junk,1)*(1:lalalala),cldpres,cldopt); shading flat; colormap(jett); colorbar; title('cld OD')
    set(gca,'ydir','reverse')  
  figure(8); clf; pcolor(ones(junk,1)*(1:lalalala),cldpres,cldde);  shading flat; colormap(jett); colorbar; title('cld DME')
    set(gca,'ydir','reverse')  

  figure(7); clf; semilogy(cldopt,cldpres); title('cld OD'); set(gca,'ydir','reverse')
  figure(8); clf; semilogy(cldde, cldpres); title('cld DME'); set(gca,'ydir','reverse')    

  %save jajaja.mat cldpres cldopt
  %pwd
  %error('ooo')
  
else
  figure(2); clf
    plot(squeeze(unique_col_frac(ibox,:,:))'); shading flat; colormap(jett); colorbar; title('cld frac')
    ylabel('vertical')
      set(gca,'ydir','reverse')
      
  figure(3); clf
    subplot(121); plot(cld_qw'); shading flat; colormap(jett); colorbar; title('cld qw')
    set(gca,'ydir','reverse')    
      ylabel('vertical')
    subplot(122); plot(cld_qi'); shading flat; colormap(jett); colorbar; title('cld qi')
    set(gca,'ydir','reverse')
    
  figure(4); clf
    subplot(121); plot(cldphase'); shading flat; colormap(jett); colorbar; title('cld phase')
    set(gca,'ydir','reverse')    
      ylabel('vertical')
    subplot(122); plot(cld_id'); shading flat; colormap(jett); colorbar; title('cld ID')
    set(gca,'ydir','reverse')
    
  figure(5); clf
    subplot(121); plot(cldopt'); shading flat; colormap(jett); colorbar; title('cld OD')
        set(gca,'ydir','reverse')
      ylabel('vertical')
    subplot(122); plot(cldde'); shading flat; colormap(jett); colorbar; title('cld eff diam')
    set(gca,'ydir','reverse')
    
  figure(6); clf; plot(cldpres'); shading flat; colormap(jett); colorbar; title('cld pres')
    set(gca,'ydir','reverse')  
end
    
disp('ret to continue to next chunk'); pause
