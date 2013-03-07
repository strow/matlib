# Matlab routines for computing allsky radiances    
# (in presence of clouds)

## SARTA : uses sarta cloudy code

Main driver file usage   
     p1 = driver_sarta_cloud_rtp(h,ha,p0,pa,*run_sarta*)     
where p0 has the ERA/ECMWF **levels** atmospheric (gas,temperature, stemp) profiles **and** cloud profiles  
and *run_sarta* is an optional arument that has the following fields   
 run_sarta.cloud = +/-1 for yes/no    
 run_sarta.klayers_code = string to klayers   
 run_sarta.sartacloud_code = string to sarta cloud executable   

#### input  
     (h,ha,p0,pa) are the usual structures read in using rtpread
|field name | Description                          |
| :-  | :- |  
|h    | usual input structure with eg channel info |  
|p0   | usual input from ECMWF or ERA levels       |  
|ha,pa| usual attributes from rtpread              |  

#### output
     p1 : output structure; if run_sarta.cloud = +1 then sarta-cloudy is run, and rcalcs used

## PCRTM : uses PCRTM cloudy code

Main driver file usage   
     p1 = driver_pcrtm_cloud_rtp(h,ha,p0,pa,*run_sarta*)     
where p0 has the ERA/ECMWF **levels** atmospheric (gas,temperature, stemp) profiles **and** cloud profiles  
and *run_sarta* is an optional arument that has the following fields   
 run_sarta.clear = +/-1 for yes/no    
 run_sarta.cloud = +/-1 for yes/no    
 run_sarta.klayers_code = string to klayers   
 run_sarta.sartaclear_code = string to sarta clear executable   
 run_sarta.sartacloud_code = string to sarta cloud executable   

#### input  
     (h,ha,p0,pa) are the usual structures read in using rtpread
|field name | Description                          |
| :-  | :- |  
|h    | usual input structure with eg channel info |  
|p0   | usual input from ECMWF or ERA levels       |  
|ha,pa| usual attributes from rtpread              |  

#### output
     p1 : output structure with cloudy calcs and extra fields containing average cloud fields etc   
Examples of extra fields     
  
|name         | size              | Description                                              |   
| :- | :- | :-: |     
| sarta_clear | [41x2700 single]  |     sarta-clr calcs                                      |    
|  rad_allsky | [41x2700 double]  |     pcrtm-cld calcs                                      |    
|  rad_clrsky | [41x2700 double]  |     pcrtm - clr                                          |    
|        ncol | [1x2700 double]   |     number of random cloud overlap columns used          |    
| sarta_clr_co2_used | [1x2700 double]   |     co2 ppmv used in sarta-clr                           |    
|     pcrtm_co2_used | [1x2700 double]   |     co2 ppmv used in pcrtm-clr                           |    
|            overlap | [1x2700 double]   |     overlap switch (3 for max. random overlap)           |    
|        pcrtm_iceOD | [1x2700 double]   |     mean ice OD   |    
|       pcrtm_iceDME | [1x2700 double]   |     mean ice DME  |    
|      pcrtm_iceCTOP | [1x2700 double]   |     mean ice CTOP |    
|      pcrtm_waterOD | [1x2700 double]   |     mean water OD |    
|     pcrtm_waterDME | [1x2700 double]   |     mean water DME |    
|    pcrtm_waterCTOP | [1x2700 double]   |     mean water CTOP |    


no columns working????????????????????????????
