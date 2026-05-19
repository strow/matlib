function [prof,orig_slabs] = driver_wrapperArbNlevs_sarta_cloud_rtp(h,ha,p,pa,run_sarta)

nlevsUNIQUE = unique(p.nlevs);
if length(nlevsUNIQUE) == 1
  [prof,orig_slabs] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
else
  prof = p;
  prof.ctype  = zeros(size(p.stemp));
  prof.ctype2 = zeros(size(p.stemp));

  prof.cngwat  = zeros(size(p.stemp));
  prof.cngwat2 = zeros(size(p.stemp));

  prof.cfrac  = zeros(size(p.stemp));
  prof.cfrac2 = zeros(size(p.stemp));
  prof.cfrac12 = zeros(size(p.stemp));

  prof.cprtop  = ones(size(p.stemp)) * -9999;
  prof.cprtop2 = ones(size(p.stemp)) * -9999;
  prof.cprbot  = ones(size(p.stemp)) * -9999;
  prof.cprbot2 = ones(size(p.stemp)) * -9999;
  prof.cpsize  = ones(size(p.stemp)) * -9999;
  prof.cpsize2 = ones(size(p.stemp)) * -9999;

  %%%%%%%%%%%%%%%%%%%%%%%%%
  if nargout == 2
    orig_slabs.ctype  = zeros(size(p.stemp));
    orig_slabs.ctype2 = zeros(size(p.stemp));
  
    orig_slabs.cngwat  = zeros(size(p.stemp));
    orig_slabs.cngwat2 = zeros(size(p.stemp));
  
    orig_slabs.cfrac  = zeros(size(p.stemp));
    orig_slabs.cfrac2 = zeros(size(p.stemp));
    orig_slabs.cfrac12 = zeros(size(p.stemp));
  
    orig_slabs.cprtop  = ones(size(p.stemp)) * -9999;
    orig_slabs.cprtop2 = ones(size(p.stemp)) * -9999;
    orig_slabs.cprbot  = ones(size(p.stemp)) * -9999;
    orig_slabs.cprbot2 = ones(size(p.stemp)) * -9999;
    orig_slabs.cpsize  = ones(size(p.stemp)) * -9999;
    orig_slabs.cpsize2 = ones(size(p.stemp)) * -9999;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%

  disp(' ')
  fprintf(1,'there are %5i different nlevs out of %5i fovs \n',length(nlevsUNIQUE),length(p.stemp));
  for ii = 1 : length(nlevsUNIQUE)
    iafovs = find(p.nlevs == nlevsUNIQUE(ii));
    fprintf(1,'  processing %5i of %5i unique nlevs (which has %5i of %5i fovs) \n',ii,length(nlevsUNIQUE),length(iafovs),length(p.stemp));
    [h,psubset] = subset_rtp_allcloudfields(h,p,[],[],iafovs);
    if nargout == 2
      orig_slabs.ctype(iafovs)  = p.ctype;
      orig_slabs.ctype2(iafovs) = p.ctype2;
    
      orig_slabs.cngwat(iafovs)  = p.cngwat;
      orig_slabs.cngwat2(iafovs) = p.cngwat2;
    
      orig_slabs.cfrac(iafovs)  = p.cfrac;
      orig_slabs.cfrac2(iafovs) = p.cfrac2;
      orig_slabs.cfrac12(iafovs) = p.cfrac12;
    
      orig_slabs.cprtop(iafovs)  = p.cprtop;
      orig_slabs.cprtop2(iafovs) = p.cprtop2;
      orig_slabs.cprbot(iafovs)  = p.cprbot;
      orig_slabs.cprbot2(iafovs) = p.cprbot2;
      orig_slabs.cpsize(iafovs)  = p.cpsize;
      orig_slabs.cpsize2(iafovs) = p.cpsize2;
  
      orig_slabs.rcalc(:,iafovs)            = p.rcalc;
      orig_slabs.sarta_rclearcalc(:,iafovs) = p.sarta_rclearcalc;
    end

    [profx,orig_slabsx] = driver_sarta_cloud_rtp(h,ha,psubset,pa,run_sarta);

    prof.ctype(iafovs)  = profx.ctype;
    prof.ctype2(iafovs) = profx.ctype2;
  
    prof.cngwat(iafovs)  = profx.cngwat;
    prof.cngwat2(iafovs) = profx.cngwat2;
  
    prof.cfrac(iafovs)  = profx.cfrac;
    prof.cfrac2(iafovs) = profx.cfrac2;
    prof.cfrac12(iafovs) = profx.cfrac12;
  
    prof.cprtop(iafovs)  = profx.cprtop;
    prof.cprtop2(iafovs) = profx.cprtop2;
    prof.cprbot(iafovs)  = profx.cprbot;
    prof.cprbot2(iafovs) = profx.cprbot2;
    prof.cpsize(iafovs)  = profx.cpsize;
    prof.cpsize2(iafovs) = profx.cpsize2;

    prof.rcalc(:,iafovs)            = profx.rcalc;
    prof.sarta_rclearcalc(:,iafovs) = profx.sarta_rclearcalc;
    
  end
end

