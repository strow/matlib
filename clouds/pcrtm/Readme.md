This directory contains routines to use the PCRTM cloud code 

Main driver file
  [p1ALL] = driver_pcrtm_cloud_rtp(h,ha,p0ALL,pa)  
where p0all has the levels atmospheric profiles and cloud profiles

input
  h     = usual input structure with eg channel info
  p0all = usual input from ECMWF or ERA levels
  ha,pa = usual attributes from rtpread

p1all = output, with cloudy calcs and extra fields containing average cloud fields etc
