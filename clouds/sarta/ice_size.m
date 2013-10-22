function dme = ice_size(T,vers)

%% function to plot ice dme as fcn of T(K), depends on version
%% code pulled out by fake_cpsize.m
%% 
%% vers 0 = Scott default
%%      1 = PCRTM 
%%      2 = Arnott
%% eg T = 200 : 280; d0 = ice_size(T,0); d1 = ice_size(T,1); d2 = ice_size(T,2); plot(T,[d0; d1; d2],'o-')

if vers == 0
  temp = T;
  icesize = [ 30,  90, 170];
  icetemp = [213, 238, 263];
  icestd  = [  4,  15,  30];
  nicesize = length(icesize);
  minicesize = 20;
  maxicesize = 200;
  indi = 1:length(T);
   ilo = indi( find(temp(indi) <= icetemp(1)) );
   cpsize(ilo) = icesize(1);
   cpstd(ilo) = icestd(1);
   ihi = indi( find(temp(indi) > icetemp(nicesize)) );
   cpsize(ihi) = icesize(nicesize);
   cpstd(ihi) = icestd(nicesize);
   for ii=1:(nicesize-1)
      jj = indi( find(temp(indi) > icetemp(ii) & ...
                      temp(indi) <= icetemp(ii+1)) );
      cpsize(jj) = icesize(ii) + (temp(jj) - icetemp(ii))*...
         (icesize(ii+1)-icesize(ii))/(icetemp(ii+1)-icetemp(ii));
      cpstd(jj) = icestd(ii) + (temp(jj) - icetemp(ii))*...
         (icestd(ii+1)-icestd(ii))/(icetemp(ii+1)-icetemp(ii));
   end
  cpsize = cpsize + cpstd;
  dme = cpsize;

elseif vers == 1
  % coefficents from S-C Ou, K-N. Liou, Atmospheric Research
  % 35(1995):127-138.
  % for computing ice cloud effective size
  c0 = 326.3;
  c1 = 12.42;
  c2 = 0.197;
  c3 = 0.0012;

  tcld = T - 273.16;

  T2 = -50;
  T2 = -60;
  boo = find(tcld < T2); 
  tcld(boo) = T2;

  T1 = -25;
  T1 = -20;
  %T1 = -00;
  boo = find(tcld > T1); 
  tcld(boo) = T1;

  dme = c0 + c1 * tcld + c2 * tcld.^2 + c3 * tcld.^3;

elseif vers == 2
  % A GCM parameterization for bimodal size spectra and ice mass removal rates in mid-latitude cirrus clouds
  % Dorothea Ivanov    David L Mitchella,    W.Patrick Arnotta,    Michael Poellotb
  % Atmospheric Research, Volumes 59–60, October–December 2001, Pages 89–113,
  % 13th International Conference on Clouds and Precipitation
  tcld = T - 273.16;
  dme = 75.3 + 0.58954 * tcld;

end
