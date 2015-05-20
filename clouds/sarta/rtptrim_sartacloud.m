function [h,ha,p,pa] = rtptrim_sartacloud(h,ha,p,pa)

%% this function is similar to Pauls rtptrim, except it works on 
%% the cloud fioleds that are produced

xfield = 'sarta_lvlZ';  %% these are heights
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfield = 'sarta_wgtI';  %% these are wgt fcns for ICE only RT
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_index_wgtpeakI';  %% these are wgt fcns peak indices for ICE
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_wgtpeakI';  %% these are wgt fcns peak indices for ICE
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_lvlDMEice';  %% these are ice DME 
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_iceOD_warn';  %% these are warnings about ice
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_lvl_iceOD_1';  %% this is where iceOD ~ 1
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'icecldX';  %% this is where we try different weighting schemes
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'icecldY';  %% this is where we try different weighting schemes
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfield = 'sarta_wgtW';  %% these are wgt fcns for WATER only RT
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_index_wgtpeakW';  %% these are wgt fcns peak indices for WATER
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_wgtpeakW';  %% these are wgt fcns peak indices for WATER
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_lvlDMEwater';  %% these are water DME 
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_waterOD_warn';  %% these are warnings about water
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'sarta_lvl_waterOD_1';  %% this is where waterOD ~ 1
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'watercldX';  %% this is where we try different weighting schemes
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

xfield = 'watercldY';  %% this is where we try different weighting schemes
if isfield(p,xfield)   
  p = rmfield(p,xfield);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDoThis = -1;
if iDoThis > 0
  xfield = 'efreq';  
  if isfield(p,xfield)   
    p = rmfield(p,xfield);
  end

  xfield = 'emis';  
  if isfield(p,xfield)   
    p = rmfield(p,xfield);
  end

  xfield = 'rho';  
  if isfield(p,xfield)   
    p = rmfield(p,xfield);
  end

  xfield = 'udef';  
  if isfield(p,xfield)   
    p = rmfield(p,xfield);
  end
end

