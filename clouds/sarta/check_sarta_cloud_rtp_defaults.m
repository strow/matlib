%% SLAB
  % run_sarta.klayers_code    = '/asl/packages/klayers/Bin/klayers_airs';
  % run_sarta.sartacloud_code = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  % run_sarta.klayers_code    = '/asl/packages/klayersV205/BinV201/klayers_airs';
  % run_sarta.sartacloud_code = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
%% SLAB

%% 100 layer
  %% see /home/sergio/klayersV205/Src_rtpV201_100layercloudamountsize
  %  run_sarta.klayers_code = ...
  %      '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';

  %% see /home/sergio/SARTA_CLOUDY/v108_Src_rtpV201_pclsam_slabcloud_hg3_100layerNEW
  %  run_sarta.sartacloud_code = ...
  %      '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';
%% 100 layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% defaults

if narginx == 4
  %% default to running sarta_cloudy
  run_sarta.clear   = -1;  %% do not run clear code
  run_sarta.cloud   = +1;  %% run cloudy code

  run_sarta.cumsum  = -1;    %% use pre-2012 cloudtop heights, without adjustments
  run_sarta.cumsum  = 9999;  %% use this in later runs eg
                             %% ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA_CLOUD_NADIR/clustbatch_make_eracloudrtp_nadir_sarta_filelist.m

  run_sarta.cfrac               = -1;  %% use random cfracs (instead of fixed fractions set by run_sarta.cfrac > 0)
  run_sarta.ice_water_separator = -1;  %% do not separate out ciwc and clwc by pressure; ie believe the NWP are correct
  run_sarta.randomCpsize        = +1;  %% keep randomizing dme for ice and water
  run_sarta.co2ppm              = 385;
  run_sarta.ForceNewSlabs       = -1;  %% keep slab clouds that are found as they are
  run_sarta.tcc                 = +1;  %% if this field exists in input structure p, then set p.cfrac = p.tcc as this is GOOOD
  
  run_sarta.waterORice = +1; % keep only water clds   %% this is for driver_sarta_cloud_rtp_onecldtest.m
  run_sarta.waterORice = -1; % keep only ice   clds   %% this is for driver_sarta_cloud_rtp_onecldtest.m
  run_sarta.waterORice = 0;  % keep both water and ice clouds ie does nothing

  run_sarta.iNew_or_Orig_CXWC2OD = -1;                %% (default) is to do OD = blah * qBlah / cc * diffZ almost PCRTM way
  
  addpath ../
  choose_klayers_sarta   %% this is for two slab only

  run_sarta.Slab_or_100layer = +1;     %% run Slab clouds
  if run_sarta.Slab_or_100layer == -1  %% which it WILL NOT be, given line above, these would be the default settings

    run_sarta.ncol    =  1;  %% number of columns for 100 layer cloud code
    run_sarta.ncol    = 25;  %% number of columns for 100 layer cloud code
    run_sarta.overlap = +3;  %% maximal random overlap for 100 layer cloud code

    %% see /home/sergio/klayersV205/Src_rtpV201_100layercloudamountsize
    run_sarta.klayers_code = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';

    %% see /home/sergio/SARTA_CLOUDY/v108_Src_rtpV201_pclsam_slabcloud_hg3_100layerNEW
    run_sarta.sartacloud_code = ...
       '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';
  end
    
elseif narginx == 5
  if ~isfield(run_sarta,'iNew_or_Orig_CXWC2OD')
    run_sarta.iNew_or_Orig_CXWC2OD = -1;  %%% stick  to OD = blah * qBlah / cc * diffZ                    Pre March 2017  DEFAULT
  end
  
  if ~isfield(run_sarta,'waterORice')
    run_sarta.waterORice = +1; % keep only water clds   %% this is for driver_sarta_cloud_rtp_onecldtest.m
    run_sarta.waterORice = -1; % keep only ice   clds   %% this is for driver_sarta_cloud_rtp_onecldtest.m
    run_sarta.waterORice = 0;  % keep both water and ice clouds ie does nothing
  end  
  if ~isfield(run_sarta,'ForceNewSlabs')
     run_sarta.ForceNewSlabs       = -1;  %% keep slab clouds that are found as they are; irrelevant for driver_sarta_cloud_rtp
  end
  if ~isfield(run_sarta,'tcc')
     run_sarta.tcc  = +1;  %% if this field exists in input structure p, then set p.cfrac = p.tcc as this is GOOOD
  end

  if ~isfield(run_sarta,'co2ppm')
    run_sarta.co2ppm = 385;
  end
  if ~isfield(run_sarta,'randomCpsize')
    run_sarta.randomCpsize = +1;
  end
  if ~isfield(run_sarta,'ice_water_separator')
    run_sarta.ice_water_separator = -1;
  end
  if ~isfield(run_sarta,'clear')
    run_sarta.clear = -1;
  end
  if ~isfield(run_sarta,'cloud')
    run_sarta.cloud = +1;
  end
  if ~isfield(run_sarta,'cumsum')
    run_sarta.cumsum = -1;    %% use pre-2012 cloudtop heights, without adjustments
    run_sarta.cumsum = 9999;  %% use this in later runs eg
                            %% ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA_CLOUD_NADIR/clustbatch_make_eracloudrtp_nadir_sarta_filelist.m
  end
  if ~isfield(run_sarta,'cfrac')
    run_sarta.cfrac = -1;  %% use random cfracs (instead of fixed fractions set by run_sarta.cfrac > 0)  
  end

  if ~isfield(run_sarta,'Slab_or_100layer')
    run_sarta.Slab_or_100layer = +1;     %% run Slab clouds
  end

  if run_sarta.Slab_or_100layer == +1    %% run 2 slab clouds
    addpath ../
    choose_klayers_sarta
  elseif run_sarta.Slab_or_100layer == -1    %% run 100layer clouds
    if ~isfield(run_sarta,'ncol')
      run_sarta.ncol =  1;  %% number of columns for 100 layer cloud code
      run_sarta.ncol = 25;  %% number of columns for 100 layer cloud code
    end
    if ~isfield(run_sarta,'overlap')
      run_sarta.overlap = +3;  %% maximal random overlap for 100 layer cloud code    
    end

    if ~isfield(run_sarta,'klayers_code')
      %% see /home/sergio/klayersV205/Src_rtpV201_100layercloudamountsize
      run_sarta.klayers_code = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';
    end

    if ~isfield(run_sarta,'sartacloud_code')
      %% see /home/sergio/SARTA_CLOUDY/v108_Src_rtpV201_pclsam_slabcloud_hg3_100layerNEW
      run_sarta.sartacloud_code = ...
       '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';
    end
    
  end             %% run_sarta.Slab_or_100layer == -1
  
end


% Min allowed cloud fraction
cmin = 0.0001;

% Max allowed cngwat[1,2]
cngwat_max = 500;

iDebugMain = +1;  %% yes debug keyboards
iDebugMain = -1;  %% no debug keyboards

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(p,'ciwc')
  error('driver_pcrtm_cloud_rtp.m requires ciwc');
elseif ~isfield(p,'clwc')
  error('driver_pcrtm_cloud_rtp.m requires clwc');
elseif ~isfield(p,'cc')
  error('driver_pcrtm_cloud_rtp.m requires cc');
elseif h.ptype ~= 0
  error('driver_pcrtm_cloud_rtp.m requires LEVELS profiles (h.ptype = 0)');
end

if run_sarta.ice_water_separator > 0
  disp('>>>>>>>> warning : setting SEPARATOR for ice and water .... initializing')
  p = convert_ice_water_separator(p,run_sarta.ice_water_separator);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
