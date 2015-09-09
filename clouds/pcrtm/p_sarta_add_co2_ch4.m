function [hxjunk,pxjunk] = p_sarta_add_co2_ch4(h,p,sarta_gas_2_6,p0ALL);

%% adds in CO2 and CH4 profiles, first to Top-Of-Current-Atmosphere (eg 1 mb for ERA) and then all the way to 0.005 mb

if h.ptype == 0
  if sum(double(h.gunit) - [21 21]') ~= 0
    error('oops gas units problem in get_sarta_clear2')
  end
  if sum(double(h.glist) - [1 3]') ~= 0
    error('oops gas list problem in get_sarta_clear2')
  end

  iDoKlayersHard = -1;
  iDoKlayersHard = +1;

  % if iDoKlayersHard == -1
  %   p_add_co2_ch4_simple    %% just adds in CO2 profiles to 0.1 mb (top of ERA and ECM atm); also add in CH4
  % end
  
  %% adds in CO2/CH4 profile, and WV and O3 and Tz upto 0.005 mb, essentially calls p_add_co2_simple
  [hxjunk,pxjunk] =  p_add_co2_ch4_complete(h,p,sarta_gas_2_6,p0ALL,iDoKlayersHard);
  
else
  disp(' OOPS cannot alter co2 profile as h.ptype is NOT 0')
end

