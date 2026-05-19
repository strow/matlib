function prof = get_sarta_cloud100layerNtimes(h,ha,prof0,pa,run_sarta)

prof = prof0;

tic

profNCOL_IP = prof;   %%%% <<<<<<<<<<<<<<<< save this, need it a lot!!!!!
h_IP          = h;

[hjunk,hajunk,pjunk,pajunk] = get_sarta_cloud100layer_klayersONLY(h,ha,prof,pa,run_sarta);

profNCOL_OP   = pjunk;
h_OP          = hjunk;
[mmjunk,nnjunk] = size(profNCOL_OP.plevs);
[unique_col_frac,ucol_num,ucol_num_same,subcol_frac] = ...
   get_subcolumn_frac_v2(length(prof.stemp), mmjunk, run_sarta.ncol, profNCOL_OP.cc',...
                                    run_sarta.overlap);

for iCol = 1 : run_sarta.ncol
  fprintf(1,'running SARTA cloud N times : subcol %3i out of %3i \n',iCol,run_sarta.ncol);
  prof = profNCOL_OP;
    
  [hX,profX] = do_subcol_cloudprofs(h_OP,prof,squeeze(subcol_frac(:,iCol,:)));

  %ijunk = 1271;  figure(6);
  %plot(prof.gas_201(:,ijunk),prof.plevs(:,ijunk),'bx-',prof.gas_202(:,ijunk),prof.plevs(:,ijunk),'rx-',...
  %     profX.gas_201(:,ijunk),profX.plevs(:,ijunk),'c',profX.gas_202(:,ijunk),profX.plevs(:,ijunk),'m')
  %title(num2str(iCol))
  %pause(0.1)

  if hX.ptype == 0
    xciwc(iCol,:,:) = profX.ciwc;
    xclwc(iCol,:,:) = profX.clwc;
  else
    xciwc(iCol,:,:) = profX.gas_201;
    xclwc(iCol,:,:) = profX.gas_202;
  end
  profRX2 = get_sarta_cloud100layer_sartaONLY(hX,ha,profX,pa,run_sarta);
  %% junkcalc(iCol,:,:) = profRX2.rcalc;    %% MEMORY HOG
  if iCol == 1
    %% slower, but more memory efficient
    [sumy,sumysqr,Nmatr] = accum_mean_std(0,0,0,profRX2.rcalc,1);
  else
    %% slower, but more memory efficient
    [sumy,sumysqr,Nmatr] = accum_mean_std(sumy,sumysqr,Nmatr,profRX2.rcalc,iCol);
  end
end

prof = profNCOL_IP;

%% prof.rcalc     = squeeze(nanmean(junkcalc,1));
%% prof.rcalc_std = squeeze(nanstd(junkcalc,1));

prof.rcalc100 = sumy./Nmatr;
junk_mean = prof.rcalc100;
%prof.rcalc_std = ...
%  real(sqrt((sumysqr - 2*junk_mean.*sumy + Nmatr.*junk_mean.*junk_mean)./(Nmatr-1)));
prof.rcalc100_std  = ...
  real(sqrt((sumysqr - 2*junk_mean.*sumy + Nmatr.*junk_mean.*junk_mean)./(Nmatr-0)));

toc
