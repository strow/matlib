function aa = cloud_mean_press(aaIN,xcumsum,icecld,watercld,plevs,ii)

aa = aaIN;

icecldX   = icecld;
watercldX = watercld;

%% these are basically the pdfs
icecldX = icecldX/nansum(icecldX);
watercldX = watercldX/nansum(watercldX);

%% these are the weights; 
icecldXW = ones(size(icecldX));%% unit weight
watercldXW = ones(size(watercldX));%% unit weight
icecldYW = ones(size(icecldX));%% unit weight
watercldYW = ones(size(watercldX));%% unit weight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% better to weight more at Top rather than Bottom - look at blahXW
%% look at cumsum from top to bottom               - look at blahYW
if plevs(1) < plevs(10)
  icecldXW = flipud(cumsum(flipud(icecldX)));    
  watercldXW = flipud(cumsum(flipud(watercldX)));    

  icecldYW = cumsum((icecldX));    
  watercldYW = cumsum((watercldX));    

else
  icecldXW = cumsum((icecldX));    
  watercldXW = cumsum((watercldX));    

  icecldYW = flipud(cumsum(flipud(icecldX)));    
  watercldYW = flipud(cumsum(flipud(watercldX)));    
end

%% this is mean cld pressure (a^0 = 1, so no weighting!!!)
aa.icecldX(ii) = nansum((plevs.*(icecldXW.^0)).*icecldX);
aa.watercldX(ii) = nansum((plevs.*(watercldXW.^0)).*watercldX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this is where we already have xcumsum of column total
%% xcumsum == run_sarta.cumsum
bonk = find(icecldYW > xcumsum,1);
if length(bonk) > 0 & xcumsum > 0 & xcumsum <= 1
  aa.icecldY(ii) = plevs(bonk);
else
  aa.icecldY(ii) = 1200;
end

bonk = find(watercldYW > xcumsum,1);
if length(bonk) > 0  & xcumsum > 0 & xcumsum <= 1
  aa.watercldY(ii) = plevs(bonk);
else
  aa.watercldY(ii) = 1200;  
end

iDebug = -1;
if iDebug > 0
  disp('setting cloud_mean_press.m')
  com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false)   
  keyboard
  [aa.icecldX(ii) aa.watercldX(ii) aa.icecldY(ii) aa.watercldY(ii)]
  figure(1); subplot(121); plot(icecldX, plevs,watercldX ,plevs); set(gca,'ydir','reverse')
  figure(1); subplot(122); plot(icecldXW,plevs,watercldXW,plevs); set(gca,'ydir','reverse')
%  figure(1); subplot(122); plot(1-icecldXW,plevs,1-watercldXW,plevs); set(gca,'ydir','reverse')
end
