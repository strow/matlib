function [head hattr prof pattr] = rtpadd_emis_DanZhou(head,hattr,prof,pattr)
%function [head hattr prof pattr] = rtpadd_emis_DanZhou(head,hattr,prof,pattr)
%
% Compute surface emissivity and add it to the RTP profile structure.
% For Land: Dan Zhou's emissivity climatology.
% For Sea : "cal_seaemis2.m" routine
%
% Necessary fields:
%   prof.satzen
%   prof.wspeed
%   prof.rtime/rlat/rlon
%   prof.landfrac
%

%  Written 17 March 2011 - Paul Schou
%  Finalised - 2013.05.08 - Breno Imbiriba


  debug = 1; %indicate that we are in debug an print out a bunch of checks

  if ~isfield(prof,'wspeed');
      error('Prof structure missing wspeed field')
  end


  % conver rtime to matlab time in order to extract month and day
  rtime_str = get_attr(pattr,'rtime');
  if isempty(rtime_str); 
    % backwards compatibility  
    rtime_str = get_attr(pattr,'L1bCM rtime'); 
  end  
  if(isempty(rtime_str))
    error('No definition of rtime in pattr');
  end
  if debug; 
    disp(['rtpadd_emis_DanZhou: rtime string = "' rtime_str '"']); 
  end

  % Usually this string is: "seconds since 0z 1 Jan YYYY"
  % Look for yyyy (four digits in the string)
  ist_year = regexp(rtime_str,'\d{4}');
  st_year = str2num(rtime_str(ist_year:ist_year+3));

  % Check if this year is AIRS or CrIS (1993/2000)
  if(st_year ~= 1999 & st_year ~= 2000)
    error(['Bad definition of rtime - not coded to handle ' num2str(st_year)]);
  end



  % Get land emissivity
  [efreq, emis] = emis_DanZhou(prof.rlat, prof.rlon, prof.rtime, st_year);

  % Get water emissivity
  [sea_nemis, sea_efreq, sea_emis] = cal_seaemis2(prof.satzen, prof.wspeed);

  % Mix them accordingly
  lgood_land = (all(emis>=0)); % good land emissivities
  lland = (prof.landfrac==1 & lgood_land); % land AND good land emis
  lsea = (prof.landfrac==0 | ~lgood_land); % ocean OR bad land emis
  lmix = ~lland & ~lsea; % the left over

  % Clean up arrays 
  prof.nemis = zeros([1, size(prof.rtime,2)]);
  prof.efreq = zeros([100,size(prof.rtime,2)]);
  prof.emis = zeros([100,size(prof.rtime,2)]);

  % Add land 
  for ifov = find(lland)
    nemis = numel(efreq);
    prof.nemis(1,ifov) = nemis; 
    prof.efreq(1:nemis,ifov) = efreq(1:nemis,1);
    prof.emis(1:nemis,ifov) = emis(1:nemis,ifov);
  end  
       
  % Add water
  for ifov = find(lsea)
    nemis = sea_nemis(1,ifov);
    prof.nemis(1,ifov) = nemis;
    prof.efreq(1:nemis,ifov) = sea_efreq(1:nemis,ifov);
    prof.emis(1:nemis,ifov) = sea_emis(1:nemis,ifov);
  end
  
  % The mixing requires attention:
  for ifov = find(lmix) 
  
    % Interpolate into land emis grid.
    nemis_sea = sea_nemis(1,ifov);
    nemis_land = numel(efreq);
    sea_emis_on_landgrid = interp1(sea_efreq(1:nemis_sea,ifov),sea_emis(1:nemis_sea,ifov), efreq(1:nemis_land,1),'linear');

    % Find the valid (non-NAN) points - use only them
    iok = find(~isnan(sea_emis_on_landgrid));
    nemis_mix = numel(iok);

    prof.nemis(1,ifov) = nemis_mix;
    prof.efreq(1:nemis_mix, ifov) = efreq(iok,1);

    % Mix both using landfrac
    lf = prof.landfrac(1,ifov);
    of = 1-lf;
    prof.emis(1:nemis_mix, ifov) = of*sea_emis_on_landgrid(iok, 1) + ...
				   lf*emis(iok,1);
  end

  % Compute Lambertian Reflectivity
  prof.nrho = prof.nemis;
  prof.rho = (1.0 - prof.emis)./3.14159265358979323846;


  % set an attribute string to let the rtp know what we have done
  pattr = set_attr(pattr,'emis',['land(' which('emis_DanZhou') ')  water(' which('cal_seaemis2') ')']);

end % Function end
