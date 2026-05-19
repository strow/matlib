function [cngwat_g_per_m2_sarta cT cB] = convert_gg_to_gm2(cT,cB,ciwc_clwc_gg_ecmwf,plevs,tlevs,spres,airslevels,airslayers,airsheights) 

%% this is the conversion from ERA/ECM total profile g/g to SARTA slab g/m2
%%
%% changes the [cT,cB,ciwc_clwc_gg_ecmwf,plevs,tlevs] from ECMWF (in g/g) to g/m2
%%
%% input
%%        cT,cB              = level number for cloud tops, cloud bottoms
%%        ciwc_clwc_gg_ecmwf = cumulative (total) cloud amt in ECMWF units (g/g)
%%        plevs,tlevs        = (could be augmented high res) pressure levels and level temps from ECMWF
%%        spres              = surface pressure from ECMWF
%% output
%%        cngwat_g_per_m2_sarta = cumulative (total) cloud amt in SARTA units (g/m2)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see /asl/packages/klayersV204/Doc/mr_from_amount.txt
%% Values of important constants: 
%%   Loschmidt = 2.6867775E+19 molecules per cm^3 (at 1 atm and 273.15 K)
%%   kAvogadro = 6.022142E+26 molecules per kilomole
%%   T0 = 273.15 K
%% In the following discussion, all profile values are layer average values.  
%% For a homogeneous air path, the relationship between column density
%% and mixing ratio is:
%%    CD_i = PP_i/kAvogadro * dz *  T0/T * Loschmidt    where
%%      MR_i is the volume mixing ratio of gas "i" expressed as the
%%        number of gas "i" molecules per total number of molecules
%%        making up the air.
%%      Ptotal is the total air pressure (atm)
%%   so PP_i = MR_i * Ptotal
%%           = the partial pressure of gas "i" (atm)
%%      CD_i is the column density of gas "i" (kilomoles/cm^2)
%%      dz is the pathlength (cm)
%%      T is the gas temperature (K)
%%
%% according to me,
%%    ---------------------------------------------------
%%    CD_i = PP_i/P0 * dz *  T0/T * Loschmidt / kAvogadro
%%    ---------------------------------------------------
%% proof 1 : by units/ PP_i/PP and T0/T have no units
%%           dz = units of cm
%%           Loschmiddt / Avogadro = molecules per cm3 / molecules per kilmole
%%                                 = kilomole per cm3
%%           Thus [] [cm] [] [mol/cm3] = kilomol/cm2
%%                                CD_i = column density of gas "i" (kilomoles/cm^2)
%% 
%% proof 2 : Po Vo = No k To  where Po, Vo, To = STP values
%%         Lo = loschmidt number = No/Vo = Po/(kTo) = 2.6867775E+19 mols/cm^3
%% also at level z, pV = N k T = n NN k T        N  = number of molecules
%%                                               NN = Avogadro number
%%                                               n  = number of moles
%%             thus  p(Adz) = n/1000 1000NA k T 
%%                   p(Adz) = n'     NN'    k T  n' = kilomoles
%%                                               NN'= kiloAvogadro
%%             thus n'/A = kilomoles/cm2 =  p dz/(NN' k T)
%% but Po Vo = No k To 
%%     p  V  = No k T           ==> 1/(kT) = No/(pV)
%%  Recall Loschmidt Lo = No/Vo ==> 1/(kT) = Lo Vo / (pV) = Lo/p Vo/V
%% Using Po Vo / To = p V / T we have  =  Vo/V = pTo/(PoT)
%%
%% Thus 1   = Lo p To
%      ---    -------
%      kT      p Po T
%%
%% which means n'     p dz     p dz  Lo p To    p   Lo  To  dz
%%             -- = ------- =  ---- --------  = -- ---- --
%%             A    NN' k T     NN' p  Po  T    Po  NN' T
%%
%% Looking at /asl/packages/klayersV204/Src/toppmv.f
%%                 PPMV = MRgg*(MDAIR/MASSF)*1E+6      MDAIR = 28.966 g/mol
%%                                                     MASSF = 18     g/mol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load airslevels.dat

global iWhichInterp  %% 0 = matlab interp1, 1 = interp1qr, set in driver_sarta_cloud_rtp.m

%% have to sanity check plevs, tlevs, spres
bad = find(plevs <= spres+100 & tlevs < 150);
good = find(plevs <= spres & tlevs > 150);
if length(bad) > 0
  tlevs0 = tlevs;
  disp('warning .. find -ve temps in the input temperature profile, temporily extrapolating them outta existence')
  tlevs(bad) = interp1(log(plevs(good)),tlevs(good),log(plevs(bad)),[],'extrap');
  %semilogy(tlevs0,plevs,'x-',tlevs,plevs,'o-'); set(gca,'ydir','reverse'); ax = axis; line([ax(1) ax(2)],[spres spres],'color','k');
  %xlim([180 320])
end

Loschmidt = 2.6867775E+19; % molecules per cm^3 (at 1 atm and 273.15 K)
kAvogadro = 6.022142E+26;  % molecules per kilomole
T0        = 273.15;        % Kelvin
P0        = 1013.5;        % mb

MDAIR = 28.966; % g/mol
MASSF = 18;     % g/mol for both water and ice!

for ii = 1 : length(cT)

  if cB(ii) > length(plevs)
    cB(ii) = length(plevs);
  end

  if cT(ii) > length(plevs)
    cT(ii) = length(plevs);
  end

  if (cT(ii) == cB(ii))
    disp('>>>>>> warning .... convert_gg_to_gm2 .... cT(ii) == cB(ii) ... quickfix else will get cngwat = NAN')
    cT(ii) = max(1,cT(ii)-2);
    cB(ii) = cT(ii)+2;
  end

  if (cT(ii) > cB(ii))
    disp('>>>>>> warning .... convert_gg_to_gm2 .... cT(ii) > cB(ii) ... quickfix else will get cngwat = INF')
    cT(ii) = max(1,cT(ii)-2);
    cB(ii) = cT(ii)+2;
  end

  %% WRONG as this is INTEGRATED amount
  %% mr = ciwc_clwc_gg_ecmwf(ii);   

  %% CORRECT  : need to go back to "mr per layer"
  mr = ciwc_clwc_gg_ecmwf(ii)/(cB(ii)-cT(ii)+1);

  if iWhichInterp == 0
    tnew = interp1(log(plevs),tlevs,log(airslevels),'spline','extrap');
  elseif iWhichInterp == 1
    tnew = interp1qr(log(plevs),tlevs,log(airslevels));
  end
  jjT = find(airslevels <= plevs(cT(ii))); jjT = min(jjT);
  jjB = find(airslevels >= plevs(cB(ii))); jjB = max(jjB);
  jj = [jjB : jjT];

  iDoSum = -1;
  if iDoSum > 0
    clear g_m2new sum_g_m2new
    g_m2new     = 0.0;
    sum_g_m2new = 0.0;
    for jjind = 1 : length(jj)
      pB = airslevels(jj(jjind));
      pT = airslevels(jj(jjind)+1);
      pnew = (pB-pT)/log(pB/pT);

      tB = tnew(jj(jjind));
      tT = tnew(jj(jjind)+1);
      slope = (tB-tT)/(log(pB)-log(pT));
      Tnew = tB - slope*(log(pB)-log(pnew));

      %hB = p2hFAST(pB,airslevels,airslayers,airsheights);
      %hT = p2hFAST(pT,airslevels,airslayers,airsheights);
      %dz = (hT - hB)*100;
      [hBT] = p2hFAST2(pB,pT,airslevels,airslayers,airsheights);
      dz = (hBT(2) - hBT(1))*100;

      ppmv = mr * (MDAIR/MASSF)*1E+6;
      pp = ppmv/1e6 * pnew;  %%% use volume mix ratio rather than ppmv
      num_kmoles_percm2 = pp/P0 * Loschmidt/kAvogadro * T0/Tnew * dz;
      %%kilomoles -> moles, cm2 -> m2
      g_m2new(jjind) = num_kmoles_percm2 * 1000 * MASSF * 10000;  
    end
    sum_g_m2new = sum(g_m2new(g_m2new >= 0));  
  end

  pB = plevs(cB(ii));
  pT = plevs(cT(ii));
  pnew = (pB-pT)/log(pB/pT);
  % pold  = pB - (pB-pT)/2;  
  % [pT pB pold pnew];

  tB = tlevs(cB(ii));
  tT = tlevs(cT(ii));
  slope = (tB-tT)/(log(pB)-log(pT));
  Tnew = tB - slope*(log(pB)-log(pnew));
  % Told  = tB - (tB-tT)/2;  
  % [tT tB Told Tnew];

  %hB = p2hFAST(pB,airslevels,airslayers,airsheights);
  %hT = p2hFAST(pT,airslevels,airslayers,airsheights);
  %dz = (hT - hB)*100;
  [hBT] = p2hFAST2(pB,pT,airslevels,airslayers,airsheights);
  dz = (hBT(2) - hBT(1))*100;

  % ppmv    = mr * (MDAIR/MASSF)*1E+6;
  % pp      = ppmv/1e6 * pold;                    %% use volume mix ratio rather than ppmv
  % num_kmoles_percm2 = pp/P0 * Loschmidt/kAvogadro * T0/Told * dz;  %% kilomole/cm2
  % g_m2old = num_kmoles_percm2*1000*MASSF*10000; %% kilomoles -> moles, cm2 -> m2

  ppmv    = mr * (MDAIR/MASSF) * 1E+6;
  pp      = ppmv/1e6 * pnew;                      %% use volume mix ratio rather than ppmv
  num_kmoles_percm2 = pp/P0 * Loschmidt/kAvogadro * T0/Tnew * dz;  %% kilomole/cm2
  g_m2new = num_kmoles_percm2*1000*MASSF*10000;   %% kilomoles -> moles, cm2 -> m2

  g_m2 = g_m2new;
  cngwat_g_per_m2_sarta(ii) = g_m2;

  %format short e
  %[mr g_m2 g_m2new sum_g_m2new]
  %format 
end

