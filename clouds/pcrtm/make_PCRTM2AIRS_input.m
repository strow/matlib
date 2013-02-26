function make_PCRTM2AIRS_input(sfcp, sfcT,newP, newT, newh2o, newo3,ncol_cld,cldnum, cld_id,cldpres, cldopt, cldde, cldphase, sfctype, sfcfreq,sfcemis,zenang, newco2, parname, usrstr, endsign)

% this code is to make input files for PCRTM over each grid-box, 
% all-sky(totally ncol_cld) + % clear-sky (totally 1)

% sfcp              surface pressure in mb
% sfcT              surface temperature in K
% newP              pressure in mb, from top to surface
% newT              temperature in k, from top to surface
% newh2o            specific humidity in g/kg, from top to surface
% newo3             ozone density in g/kg, from top to surface
% ncol_cld          number of sub-columns
% cldnum            number of cloud layers of each subcolumn
% cldpre            cloud top pressure in mb of each cloud layer
% cldopt            cloud visible optical depth of each cloud layer
% cldde             cloud effective radius of each cloud layer
% cldphase          cloud phase of each cloud layer
% sfctype           surface type defined in IGBP
% sfcfreq           frequencies at which user give the surface emissivity
% sfcemis           user given surface emissivity
% if sfctype <=0    use emissivity defined by sfcemis
% zenang            viewing zenith angle
% newco2            new co2 volume mixing ratio for scaling co2 profile in ppmv
% parname           file name for PCRTM input file
% usrstr            a string
% endsign           a flag to indicate whether is the end of input or not

% error checking
if nargin ~= 21
	disp('wrong number of inputs');
	return;
end


nang =length(zenang);

wid = fopen(parname, 'a');

% here we set ncol_cld +1, because we add a clear-sky input
fprintf (wid, '%s  %5d %5d %5.1f %6.1f %5d %10.3f\n', usrstr, ncol_cld +1, nang, sfcT, sfcp, sfctype,newco2);

for ilev =1:length(newP)          
         
     fprintf (wid, '%5d %15f %15f %15e %15e \n', ilev, newP(ilev), newT(ilev), newh2o(ilev), newo3(ilev));
end


if sfctype<=0
    idx = find(isnan(sfcemis)==0);  

   if ~isempty(idx)
        fprintf (wid, '%d \n', length(idx));
        for iw = 1:length(idx)
             fprintf (wid, '%10.3f %8.3f\n', sfcfreq(idx(iw)), sfcemis(idx(iw)));
        end
   else
       sfctype =17; % simply set default surface type as ocean
    end
end


% viewing angles
for i = 1:nang
     fprintf (wid, '%8.3f \n', zenang(i));
end

% all-sky cloud parameters  
 
for icol = 1:ncol_cld
    
    fprintf (wid,'%d \n', cldnum(icol));
    
    for  i=1:cldnum(icol)
      
        fprintf (wid, '%3d %8.2f %8.4f %8.2f %2d \n', cld_id(icol,i), cldpres (icol,i), cldopt(icol,i), cldde(icol,i), cldphase(icol,i));
    end
end

% for clear-sky, no cloud parameters
fprintf (wid,'%d \n', 0);

if (endsign ==0 )

   fprintf(wid, '%s\n', '    1');
   
else
   fprintf(wid, '%s\n', '    0');
end

fclose(wid);
