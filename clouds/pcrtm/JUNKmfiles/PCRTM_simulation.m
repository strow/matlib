

clear all;

% PCRTM executable path
ppath = '/raid_backup/xiuchen/for_Sergio_AIRS/PCRTM_V2.1_for_AIRS/code_changed/Run/';

% output path for PCRTM output
outpath = '/raid01/xiuchen/pcrtm_test/';

% number of subcolumn
ncol =50;

% overlap =1 maximum overlap
% overlap =2 random overlap
% overlap =3 maximum-random overlap
overlap =3;

               
load xianglei_2012_05_01_00hrs     

lev = p.plevs(:,1);
nlev = length(lev);

% convert from kg/kg to g/kg
q  = p.gas_1 * 1e3;
nboxes = size(q,2);

% convert from kg/kg to g/kg
o3 = p.gas_3 * 1e3;
TT = p.ptemp;

cc = p.cc;
clwc = p.clwc;
ciwc = p.ciwc;

Ps = p.spres;
Ts = p.stemp;

% set default co2 as 385.85 ppmv, it can be replaced with new co2.
co2 = zeros(1,nboxes) + 385.85;
%co2 = p.co2ppmv;

% the given surface spectral emissivity will be interpolated onto 425
% mon-frequencies in PCRTM
efreq = p.efreq;
emis  = p.emis;

% convert to degree
zen_ang = acosd(p.satzen);
 

% set surface type as 0, so PCRTM will use user defined surface spectral emissivity  
sfc = zeros(1,nboxes);


parname = [outpath,'xianglei_2012_05_01_00hrs.par'];

unix(['rm -f ', parname]);
                   
if nboxes >0   
        
        for ibox =1:nboxes
             % set cloud fraction to 0 for cloud below surface pressure
             clear idp ;
             idp = find (lev> Ps(ibox));
             cc(idp,ibox) = 0.0;
        end
                 
        [rad_allsky rad_clrsky] = PCRTM_compute_for_AIRS_spectra(nboxes, nlev, ncol, overlap, lev, clwc, ciwc, cc, TT, q, o3, Ps, Ts, sfc,efreq,emis,zen_ang,co2,parname,ppath);
        % rad_allsky is in W per m^2 per sr per cm^-1
        % rad_clrsky is in W per m^2 per sr per cm^-1
        
        
        unix(['rm -f  ',parname, '.par.out']);                       
        unix(['rm -f  ',parname, '.par']);  
        
       % vid = fopen([parname,'_spec.ieee'],'ab','l');  
       % for ibox =1:nboxes
             
       %      fwrite (vid, rad_allsky(:,ibox), 'float');
       %      fwrite (vid, rad_clrsky(:,ibox), 'float');
       % end
   
       % fclose(vid);        
        

end
                 
           