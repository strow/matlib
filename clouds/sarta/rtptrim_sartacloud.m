function [h,ha,p,pa] = rtptrim_sartacloud(h,ha,p,pa)

%% this function is similar to Pauls rtptrim, except it works on 
%% the cloud fioleds that are produced

iF = 0;
clear list

list = {'sarta_lvlZ'};  %% heights
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ice
list = {'sarta_wgtI','sarta_index_wgtpeakI','sarta_wgtpeakI','sarta_lvlDMEice',...
        'sarta_iceOD_warn','sarta_lvl_iceOD_1','icecldX','icecldY'};
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% water
%% ice
list = {'sarta_wgtW','sarta_index_wgtpeakW','sarta_wgtpeakW','sarta_lvlDMEwater',...
        'sarta_waterOD_warn','sarta_lvl_waterOD_1','watercldX','watercldY'};
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p,xfield)
    iF = iF + 1;
    p = rmfield(p,xfield);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDoThis = -1;
if iDoThis > 0
  list = {'efreq','emis','rho','udef'};
  for ii = 1 : length(list)
    xfield = list{ii};
    if isfield(p,xfield)
      iF = iF + 1;
      p = rmfield(p,xfield);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

