This directory contains routines to use the PCRTM cloud code 

Main driver file
  [p1ALL] = driver_pcrtm_cloud_rtp(h,ha,p0ALL,pa)  
where p0all has the levels atmospheric profiles and cloud profiles

input
  h     = usual input structure with eg channel info
  p0all = usual input from ECMWF or ERA levels
  ha,pa = usual attributes from rtpread

p1all = output, with cloudy calcs and extra fields containing average cloud fields etc
           sarta_clear: [41x2700 single]       %% sarta - clr
            rad_allsky: [41x2700 double]       %% pcrtm - cld
            rad_clrsky: [41x2700 double]       %% pcrtm - clr
                  ncol: [1x2700 double]        
    sarta_clr_co2_used: [1x2700 double]
        pcrtm_co2_used: [1x2700 double]
               overlap: [1x2700 double]        %% these are mean cloud profiles used in the ncol simulations
           pcrtm_iceOD: [1x2700 double]
          pcrtm_iceDME: [1x2700 double]
         pcrtm_iceCTOP: [1x2700 double]
         pcrtm_waterOD: [1x2700 double]
        pcrtm_waterDME: [1x2700 double]
       pcrtm_waterCTOP: [1x2700 double]
