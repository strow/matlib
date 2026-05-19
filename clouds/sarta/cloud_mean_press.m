function aa = cloud_mean_press(aaIN,xcumsum,icecld,watercld,plevs,ii)

%% input
%%   xcumsum = run_sarta.cumsum so can be -9999,-1,-(0-1),+(0-1),1-9998,9999
%% output
%%   aa.icecldX,aa.watercldX = mean(CIWC), mean(CLWC) pressure level
%%   aa.icecldY,aa.watercldY = pressure level where normalized CIWC/CLWC exceed xcumsum if 0 < xcumsum < 1
%%                             else set to 1200 mb

aa = aaIN;

icecldX   = icecld;
watercldX = watercld;

%% these are basically the pdfs of the CIWC and CLWC profiles ....
icecldX = icecldX/nansum(icecldX);
watercldX = watercldX/nansum(watercldX);

%% these are the weights; 
icecldXW = ones(size(icecldX));     %% unit weight
watercldXW = ones(size(watercldX)); %% unit weight
icecldYW = ones(size(icecldX));     %% unit weight
watercldYW = ones(size(watercldX)); %% unit weight

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
aa.icecldX(ii)   = nansum((plevs.*(icecldXW.^0)).*icecldX);
aa.watercldX(ii) = nansum((plevs.*(watercldXW.^0)).*watercldX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if xcumsum <= 0 & xcumsum > -1
  xcumsum = rand(size(xcumsum));
end

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
  %com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false)   
  %keyboard
  [aa.icecldX(ii) aa.watercldX(ii) aa.icecldY(ii) aa.watercldY(ii)]
  figure(1); subplot(121); plot(icecldX, plevs,'b',watercldX ,plevs,'r'); set(gca,'ydir','reverse')
  figure(1); subplot(122); plot(icecldXW,plevs,'b',watercldXW,plevs,'r'); set(gca,'ydir','reverse')
%  figure(1); subplot(122); plot(1-icecldXW,plevs,1-watercldXW,plevs); set(gca,'ydir','reverse')
  error('in cloud_mean_press.m')
end
