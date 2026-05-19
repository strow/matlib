function [watercld,icecld,plevs,aa,ptemp,cut440,p440] = new_style_smooth_cc_ciwc_clwc_to_water_ice_profile(xcumsum,profX,ice_water_separator,aa,ii);

%%% old_style_smooth_cc_ciwc_clwc_to_water_ice_profile been superseded by this new code

global iWhichInterp  %% 0 = matlab interp1, 1 = interp1qr, set in driver_sarta_cloud_rtp.m

pSEPARATE = 440;
if ice_water_separator == -1
  %% default ISCCP
  pSEPARATE = 440;
elseif ice_water_separator == 0
  %% no separation at all
  pSEPARATE = 9999;
elseif ice_water_separator > 100
  pSEPARATE = ice_water_separator;

elseif ice_water_separator == 1 | ice_water_separator == 2
  X = [-60 0 +60]; Y = [6 9 6]; P = polyfit(X,Y,2);
  X1 = profX.rlat;
  Y1 = polyval(P,X1); 
  pSEPARATEX = Y1;
  plevs = profX.plevs(:,ii);
  if ~isfield(profX,'palts')
    palts = p2h(plevs);
  else
    palts = profX.palts(:,ii);
  end
  boo = find(isfinite(plevs) & isfinite(palts));
  if iWhichInterp == 0
    pSEPARATE = interp1(palts(boo),log(plevs(boo)),1000*Y1(ii));
  elseif iWhichInterp == 1
    pSEPARATE = interp1qr(palts(boo),log(plevs(boo)),1000*Y1(ii));
  end
  pSEPARATE = exp(pSEPARATE);
  if pSEPARATE < 440
    pSEPARATE = 440;
  end

  ciwc  = profX.ciwc(:,ii);
  plevs = profX.plevs(:,ii);
  plot(ciwc,plevs);
  if sum(ciwc) > eps
    nonzero_ciwc = min(find(ciwc > 0));
    if nonzero_ciwc < length(plevs)
      nonzero_ciwc = nonzero_ciwc + 1;
    end
    nonzero_ciwc = plevs(nonzero_ciwc);
    pSEPARATE = max(nonzero_ciwc,pSEPARATE);
  end
  
end

if pSEPARATE < 440
  pSEPARATE = 440;   %% else the cut off at the tropics is too high
end
p440 = pSEPARATE;

%% slabs can be resolved better if there are more points
plevs = profX.plevs(:,ii);
if length(plevs) < 80
  plevsX = (plevs(1:end-1) + plevs(2:end))/2;
  plevs = sort([plevs; plevsX]);
end
if length(plevs) < 80
  plevsX = (plevs(1:end-1) + plevs(2:end))/2;
  plevs = sort([plevs; plevsX]);
end

if plevs(1) < plevs(10)
  cut440 = find(plevs >= pSEPARATE,1);
else
  cut440 = find(plevs <= pSEPARATE,1);
end

if iWhichInterp == 0
  watercld = interp1(log10(profX.plevs(:,ii)),profX.clwc(:,ii),log10(plevs));
  icecld   = interp1(log10(profX.plevs(:,ii)),profX.ciwc(:,ii),log10(plevs));
  ptemp    = interp1(log10(profX.plevs(:,ii)),profX.ptemp(:,ii),log10(plevs));
elseif iWhichInterp == 1
  watercld = interp1qr(log10(profX.plevs(:,ii)),profX.clwc(:,ii),log10(plevs));
  icecld   = interp1qr(log10(profX.plevs(:,ii)),profX.ciwc(:,ii),log10(plevs));
  ptemp    = interp1qr(log10(profX.plevs(:,ii)),profX.ptemp(:,ii),log10(plevs));
end

%figure(2); clf;
%  subplot(121); semilogy(icecld,plevs,'b',watercld,plevs,'r'); set(gca,'ydir','reverse');
%    ax = axis; axis([ax(1) ax(2) 10 1000])
%  subplot(122); semilogy(profX.ptemp(:,ii),profX.plevs(:,ii),'bx-',ptemp,plevs,'r'); set(gca,'ydir','reverse');  
%    ax = axis; axis([220 260 10 1000])
    
if ~exist('aa','var')
  %aa = [];
  aa = struct;  
end
aa = cloud_mean_press(aa,xcumsum,icecld,watercld,plevs,ii);

