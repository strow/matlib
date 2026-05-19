function p = convert_ice_water_separator(p,pSEPARATE)

p0 = p;

%% this function mimics PCRTM in that it converts everything ABOVE psperate to ice, and below to water

iMethod = 1; %% this is simplest
iMethod = 2; %% this is harder, uses PCRTM ideassimplest

R = 287.05;
g = 9.806;

nboxes = length(p.stemp);

if iMethod == 1
  %% this is simplest
  for ibox = 1 : nboxes
    xcc  = p.cc(:,ibox);
    xice = p.ciwc(:,ibox);
    xwater = p.clwc(:,ibox);
    xpress = p.plevs(:,ibox);

    oo = find(xpress <= pSEPARATE);
      p.ciwc(oo,ibox) = xice(oo) + xwater(oo);
      p.clwc(oo,ibox) = 0;
    oo = find(xpress > pSEPARATE);
      p.clwc(oo,ibox) = xice(oo) + xwater(oo);
      p.ciwc(oo,ibox) = 0;

  end

elseif iMethod == 2
  %% this is complicated, uses PCRTM ideas

  % coefficents from S-C Ou, K-N. Liou, Atmospheric Research
  % 35(1995):127-138.
  % for computing ice cloud effective size
  c0 = 326.3;
  c1 = 12.42;
  c2 = 0.197;
  c3 = 0.0012;

  for ibox = 1 : nboxes
    xcc  = p.cc(:,ibox);
    xice = p.ciwc(:,ibox);
    xwater = p.clwc(:,ibox);
    Px = p.plevs(:,ibox);
    TT = p.ptemp(:,ibox);
    tcld = TT - 273.15;
      tcld(tcld < -50) = -50;
      tcld(tcld > -25) = -25;
    cldde_water = 20 * ones(size(Px));
    cldde_ice   = c0 + c1 * tcld + c2 * tcld.^2 + c3 * tcld.^3;

    Ttmp = 0.5*(TT(1:end-1) + TT(2:end));
    scaleH = R*Ttmp/g/1000;
    dz = scaleH .* log(Px(2:end)./Px(1:end-1));
    Z = zeros(size(Px));
    Z(2:end) = Z(2:end) + cumsum(dz);

    laZ = [abs(diff(Z)); 0];

    %%%%%%%%%%%%%%%%%%%%%%%%%
    qi = xice ./TT .* Px * 100/R*1e3; %change ice water content from kg/kg to g/m^3
    ice_opt = (0.003448 + 2.431./cldde_ice) .*qi ./ xcc .* laZ *1e3;

    qw = xwater ./ TT .*Px *100/R *1e3;  %change liquid water content from kg/kg to g/m^3
    water_opt = 3 * qw ./ cldde_water ./ xcc .* laZ *1e3;

    total_opt(:,ibox) = ice_opt + water_opt;

    oo = find(Px > pSEPARATE);
      %% make all this water
      qwx = zeros(size(xcc));
      qwx(oo) = total_opt(oo,ibox) .* cldde_water(oo) .* xcc(oo) /(3 * 1e3) ./laZ(oo);
      p.clwc(oo,ibox) = qwx(oo) .* TT(oo) * R ./ Px(oo) /1e3/100;
      p.ciwc(oo,ibox) = 0;

    oo = find(Px <= pSEPARATE);
      %% make all this ice
      qix = zeros(size(xcc));
      qix(oo) = (total_opt(oo,ibox) - 0.003448) .* cldde_ice(oo) .* xcc(oo) /(2.431 * 1e3) ./laZ(oo);
      p.ciwc(oo,ibox) = qix(oo) .* TT(oo) * R ./  Px(oo) /1e3/100;
      p.clwc(oo,ibox) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% >>>>>  this is just a sanity check, that the ODs are about the same!!!!!!! <<<<<<<<<<<<<<
    xice = p.ciwc(:,ibox);
    xwater = p.clwc(:,ibox);

    qi = xice ./TT .* Px * 100/R*1e3; %change ice water content from kg/kg to g/m^3
    ice_opt = (0.003448 + 2.431./cldde_ice) .*qi ./ xcc .* laZ *1e3;

    qw = xwater ./ TT .*Px *100/R *1e3;  %change liquid water content from kg/kg to g/m^3
    water_opt = 3 * qw ./ cldde_water ./ xcc .* laZ *1e3;

    total_opt_new(:,ibox) = ice_opt + water_opt;

  end    %% loop over ibox

boo = find(p.ciwc < 0 | isnan(p.ciwc)); p.ciwc(boo) = 0;
boo = find(p.clwc < 0 | isnan(p.clwc)); p.clwc(boo) = 0;

end
