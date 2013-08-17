function prof = set_fracs_deffs(head,pin,profX,cmin,cngwat_max,cfracSet);

prof = pin;

% Replace cfrac info

nprof = length(prof.stemp);

[cfracw, cfraci] = total_cfrac(profX.plevs,profX.cc,profX.clwc,profX.ciwc);
tcc = profX.cfrac;

[prof.cfrac, prof.cfrac2, prof.cfrac12] = fake_cfracs(tcc, cfracw, cfraci, prof.ctype, prof.ctype2);

prof.clwc = profX.clwc;
prof.ciwc = profX.ciwc;
prof.cc = profX.cc;

% Compute cloud temperature
pavg = 0.5*(prof.cprtop + prof.cprbot);
ibad = find(prof.cfrac == 0);
pavg(ibad) = 500; % safe dummy value
tavg1 = rtp_t_at_p(pavg, head, prof);
pavg = 0.5*(prof.cprtop2 + prof.cprbot2);
ibad = find(prof.cfrac2 == 0);
pavg(ibad) = 500; % safe dummy value
tavg2 = rtp_t_at_p(pavg, head, prof);
clear ibad

% Replace cpsize and cpsize2
iceflag = zeros(1,nprof);
ii = find(prof.ctype == 201);
iceflag(ii) = 1;
prof.cpsize  = fake_cpsize(tavg1, iceflag, 1);
%

iceflag = zeros(1,nprof);
ii = find(prof.ctype2 == 201);
iceflag(ii) = 1;
prof.cpsize2 = fake_cpsize(tavg2, iceflag, 1);
clear iceflag tavg1 tavg2

% Remove cloud fractions less than some minimum
hcmin = 0.5*cmin;

%%%%%%%%%%%%%%%%%%%%%%%%%
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

