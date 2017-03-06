addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER

dirMERRA = '/home/sbuczko1/testoutput/MERRA-2/2011/03/';
fname_rad_2d = [dirMERRA 'MERRA2_400.tavg1_2d_rad_Nx.20110311.nc4'];
fname_asm_3d = [dirMERRA 'MERRA2_400.tavg3_3d_asm_Nv.20110311.nc4'];
fname_flx_2d = [dirMERRA 'MERRA2_400.tavg1_2d_flx_Nx.20110311.nc4'];
fname_slv_2d = [dirMERRA 'MERRA2_400.tavg1_2d_slv_Nx.20110311.nc4'];
fname_cld_3d = [dirMERRA 'MERRA2_400.tavg3_3d_cld_Nv.20110311.nc4'];

s_rad_2d = read_netcdf_lls(fname_rad_2d);
s_asm_3d = read_netcdf_lls(fname_asm_3d);
s_flx_2d = read_netcdf_lls(fname_flx_2d);
s_slv_2d = read_netcdf_lls(fname_slv_2d);
s_cld_3d = read_netcdf_lls(fname_cld_3d);