function prof = set_fracs_deffs(head,pin,profX,cmin,cngwat_max,cfracSet,randomCpsize);

%% before Aug 2013 randomCpsize was effectively 1
%% randomCpsize = 1    ==> random dme_water (centered about 20 um), dme_ice = based on Scott Tcld, randomized
%%                20   ==> dme_water FIXED at 20 um, dme_ice = based on KNLiou Tcld, randomized
%%                -1   ==> dme_water based on MODIS climatology, dme_ice = based on Scott Tcld, randomized
%%              9999   ==> dme water based on MODIS climatology, dme_ice = based on KNLiou Tcld

prof = pin;

% Replace cfrac info

nprof = length(prof.stemp);

[cfracw, cfraci] = total_cfrac(profX.plevs,profX.cc,profX.clwc,profX.ciwc);

oo = find(profX.cfrac < 0); 
if length(oo) > 0 
  profX.cfrac(oo) = 0; 
  fprintf(1,'set_fracs_deffs : %5i profX.cfrac < 0 \n',length(oo))
end
oo = find(profX.cfrac > 1); 
if length(oo) > 0 
  profX.cfrac(oo) = 1;
  fprintf(1,'set_fracs_deffs : %5i profX.cfrac > 1 \n',length(oo))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% tcc = profX.cfrac;     %% bad, before Aug 2015
%%%%%% tcc = profX.cfrac;     %% bad, before Aug 2015
%%%%%% tcc = profX.cfrac;     %% bad, before Aug 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(profX,'tcc')
  error('looking for tcc from ERA/ECMWF/MERRA')
end
tcc = profX.tcc;              %% should be much better
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oo = find(cfracw < 0); 
if length(oo) > 0 
  cfracw(oo) = 0;
  fprintf(1,'set_fracs_deffs : %5i cfracw < 0 \n',length(oo))
end
oo = find(cfracw > 1); 
if length(oo) > 0 
  cfracw(oo) = 1;
  fprintf(1,'set_fracs_deffs : %5i cfracw > 1 \n',length(oo))
end
oo = find(cfraci < 0); 
if length(oo) > 0 
  cfraci(oo) = 0;
  fprintf(1,'set_fracs_deffs : %5i cfraci < 0 \n',length(oo))
end
oo = find(cfraci > 1); 
if length(oo) > 0 
  cfraci(oo) = 1;
  fprintf(1,'set_fracs_deffs : %5i cfraci > 1 \n',length(oo))
end

[prof.cfrac, prof.cfrac2, prof.cfrac12] = fake_cfracs(tcc, cfracw, cfraci, prof.ctype, prof.ctype2);

prof.clwc = profX.clwc;
prof.ciwc = profX.ciwc;
prof.cc = profX.cc;

prof = set_cpsize(head, prof,randomCpsize); %%% this looks new!! but really is script code that has now been put into a subroutine

% Remove cloud fractions less than some minimum
hcmin = 0.5*cmin;

if cfracSet > 0
  iPos = find(prof.cngwat > 0 & prof.cfrac > 0);
  %plot(prof.cfrac(iPos)); disp('ret'); pause
  prof.cfrac(iPos) = cfracSet;
  iPos = find(prof.cngwat2 > 0 & prof.cfrac2 > 0);
  %plot(prof.cfrac2(iPos)); disp('ret'); pause
  prof.cfrac2(iPos) = cfracSet;
end

ix=find(prof.cfrac < cmin);
prof.cfrac(ix)  = 0;
prof.cfrac12(ix)= 0;
prof.cngwat(ix) = 0;
prof.ctype(ix)  = -1;
prof.cprtop(ix) = -9999;
prof.cprbot(ix) = -9999;

ix=find(prof.cfrac2 < cmin);
prof.cfrac2(ix)  = 0;
prof.cfrac12(ix) = 0;
prof.cngwat2(ix) = 0;
prof.ctype2(ix)  = -1;
prof.cprtop2(ix) = -9999;
prof.cprbot2(ix) = -9999;

ix = find(prof.cfrac12 >= hcmin & prof.cfrac12 < cmin);
prof.cfrac12(ix) = cmin;
ix = find(prof.cfrac12 < hcmin);
prof.cfrac12(ix) = 0;
junk = prof.cfrac(ix) + prof.cfrac2(ix);
ii = ix( find(junk > 1) );
ii1 = ii( find(prof.cfrac(ii) > prof.cfrac2(ii)) );
ii2 = setdiff(ii,ii1);
prof.cfrac(ii1) = prof.cfrac(ii1)-hcmin;
prof.cfrac2(ii2) = prof.cfrac2(ii2)-hcmin;

ix = find(prof.cngwat > cngwat_max);
prof.cngwat(ix) = cngwat_max;
ix = find(prof.cngwat2 > cngwat_max);
prof.cngwat2(ix) = cngwat_max;

