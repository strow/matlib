# Matlab routines for computing allsky radiances    
# (in presence of clouds)

## SARTA : uses sarta cloudy code

## PCRTM : uses PCRTM cloudy code

Main driver file usage
  p1ALL = driver_pcrtm_cloud_rtp(h,ha,p0,pa)     
where p0 has the ERA/ECMWF **levels** atmospheric (gas,temperature, stemp) profiles **and** cloud profiles  

#### input
  |-----|--------------------------------------------|
  |h    | usual input structure with eg channel info |
  |p0   | usual input from ECMWF or ERA levels       |
  |ha,pa| usual attributes from rtpread              |

#### output
          |p1|  output |     with cloudy calcs and extra fields containing average cloud fields etc  |
          |  |         |     examples shown below                                                    |
          | sarta_clear | [41x2700 single]       %% sarta - clr | |
          |  rad_allsky| [41x2700 double]       %% pcrtm - cld  | |
          |  rad_clrsky| [41x2700 double]       %% pcrtm - clr  | |
          |        ncol| [1x2700 double]        | |
   | sarta_clr_co2_used| [1x2700 double]        | |
   |     pcrtm_co2_used| [1x2700 double]        | |
   |            overlap| [1x2700 double]        %% these are mean cloud profiles used in the ncol simulations | |
   |        pcrtm_iceOD| [1x2700 double]        |  |
   |       pcrtm_iceDME| [1x2700 double]        |  |
   |      pcrtm_iceCTOP| [1x2700 double]        |  |
   |      pcrtm_waterOD| [1x2700 double]        |  |
   |     pcrtm_waterDME| [1x2700 double]        |  |
   |    pcrtm_waterCTOP| [1x2700 double]        |  |

