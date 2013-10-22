function dme = water_size(T,vers)

%% function to plot water dme as fcn of T(K), depends on version
%% code pulled out by fake_cpsize.m
%% 
%% vers 0 = Scott default
%%      1 = PCRTM 
%% eg T = 200 : 280; d0 = water_size(T,0); d1 = water_size(T,1); d2 = water_size(T,2); plot(T,[d0; d1; d2],'o-')

if vers == 0
  temp = T;
  watsize = [ 15,   18,  21];
  wattemp = [213,  243, 273];
  watstd  = [  1,  1.5,   2];
  nwatsize = length(watsize);
  minwatsize = 13;
  maxwatsize = 26;

  indw = 1:length(T);

   ilo = indw( find(temp(indw) <= wattemp(1)) );
   cpsize(ilo) = watsize(1);
   cpstd(ilo) = watstd(1);
   ihi = indw( find(temp(indw) > wattemp(nwatsize)) );
   cpsize(ihi) = watsize(nwatsize);
   cpstd(ihi) = watstd(nwatsize);
   for ii=1:(nwatsize-1)
      jj = indw( find(temp(indw) > wattemp(ii) & ...
                      temp(indw) <= wattemp(ii+1)) );
      cpsize(jj) = watsize(ii) + (temp(jj) - wattemp(ii))*...
         (watsize(ii+1)-watsize(ii))/(wattemp(ii+1)-wattemp(ii));
      cpstd(jj) = watstd(ii) + (temp(jj) - wattemp(ii))*...
         (watstd(ii+1)-watstd(ii))/(wattemp(ii+1)-wattemp(ii));
   end
  cpsize = cpsize + cpstd;
  dme = cpsize;

elseif vers == 1
  dme = 20 * ones(size(T));

end
