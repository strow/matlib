function [ciwc_clwc_gg_ecmwf] = convert_gm2_to_gg(pT, pB, cngwat_g_per_m2, plevs, tlevs, nlevs, aLH) 

%% this is the conversion from SARTA slab g/m2 to ERA/ECM total profile g/g 
%%
%% changes the [cT,cB,ciwc_clwc_gg_ecmwf,plevs,tlevs] from ECMWF (in g/g) to g/m2
%%
%% input
%%        pT,pB            = pressure for cloud tops, cloud bottoms
%%        cngwat_g_per_m2  = slab (total) cloud amt in SARTA units (g/m2)
%%        plevs,tlevs      = pressure levels and level temps from ECMWF
%%        aLH              = structure containing aLH.airslevels, aLH.airsheights
%% output
%%        ciwc_clwc_gg_ecmwf = cumulative (total) cloud amt in ECMWF units (g/g)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% if we multiply by mass per kilomole of gas, we get g/m2
%%
%% Looking at /asl/packages/klayersV204/Src/toppmv.f
%%                 PPMV = MRgg*(MDAIR/MASSF)*1E+6      MDAIR = 28.966 g/mol
%%                                                     MASSF = 18     g/mol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airslevels  = aLH.airslevels;
airsheights = aLH.airsheights;

Loschmidt = 2.6867775E+19; % molecules per cm^3 (at 1 atm and 273.15 K)
kAvogadro = 6.022142E+26;  % molecules per kilomole
T0        = 273.15;        % Kelvin
P0        = 1013.5;        % mb

MDAIR = 28.966; % g/mol
MASSF = 18;     % g/mol for both water and ice!

for ii = 1 : length(pT)

  cT = find(plevs >= pT,1);
  cB = find(plevs >= pB,1);

  if length(cT) == 0
    cT = 1;
  end
  if length(cB) == 0
    cB = nlevs;
  end

  if cB(ii) > length(plevs)
    cB(ii) = length(plevs);
  end

  if cT(ii) > length(plevs)
    cT(ii) = length(plevs);
  end

  if (cT(ii) == cB(ii))
    disp('>>>>>> warning .... convert_gm2_to_gg .... cT(ii) == cB(ii) ... quickfix else will get cngwat = NAN')
    cT(ii) = max(1,cT(ii)-2);
    cB(ii) = cT(ii)+2;
  end

  if (cT(ii) > cB(ii))
    disp('>>>>>> warning .... convert_gm2_to_gg .... cT(ii) > cB(ii) ... quickfix else will get cngwat = INF')
    cT(ii) = max(1,cT(ii)-2);
    cB(ii) = cT(ii)+2;
  end

  tnew = interp1(log(plevs),tlevs,log(airslevels),'spline','extrap');
  jjT = find(airslevels <= plevs(cT(ii))); jjT = min(jjT);
  jjB = find(airslevels >= plevs(cB(ii))); jjB = max(jjB);
  jj = [jjB : jjT];

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

  hB = p2hFAST(pB,airslevels,airsheights);
  hT = p2hFAST(pT,airslevels,airsheights);
  dz = (hT - hB)*100;

  g_m2new           = cngwat_g_per_m2(ii);
  num_kmoles_percm2 = g_m2new/(1000*MASSF*10000); %% kilomole/cm2
  pp                = num_kmoles_percm2 * P0 * kAvogadro/Loschmidt * Tnew/T0 / dz;  %% partial pressure of cloud
  ppmv              = pp*1e6 / pnew;              %% ppmv
  mr                = ppmv/1e6 * (MASSF/MDAIR);   %% mix ratio at each point

  ciwc_clwc_gg_ecmwf(ii) = mr;
  
end
