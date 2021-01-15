% NAME
%   aslmap - global map for equal area trapezoids
%
% SYNOPSIS
%   fh = equal_area_map(fn, glat, glon, gval, latlim, lonlim, opts)
%
% INPUTS
%   fn    - matlab figure number
%   glat  - n+1 vector of latitude boundaries
%   glon  - m+1 vector of longitude boundaries
%   gval  - n x m array of map data values
%   latlim - 1 x 2 vector of maplatlimit  (display only)
%   lonlim - 1 x 2 vector of maplonlimit  (display only)
%   opts: optional structure with more map properties, possible values are:
%     opts.color =  'k'; for black lat/lon lines, etc., or 3 element vector [r g b]
%     opts.title = 'A title';
%     opts.caxis = [minvalue maxvalue];  Caxis limits
%     opts.cmap  =  jet; Some colormap, can be default  (ie. opt.cmap = 'default')
%     opts.titlesize = 14; Size of title fonts
%
% OUTPUT
%   fh    - figure handle
%
% DISCUSSION
%
function fh = aslmap(fn, glat, glon, gval, latlim, lonlim, varargin);
   
   if  size(varargin) == [1 1]
      do_opts = 1;
      opts = varargin{1};
   else
      do_opts = 0;
   end
   
% NaN grid extension
   [m, n] = size(gval);
   gtmp = NaN(m+1,n+1);
   gtmp(1:m, 1:n) = gval;

% 2d lat and lon arrays
   glat = glat(:) * ones(1, n+1);
   glon = ones(m+1,1) * glon(:)';

   fh = figure(fn);  clf
   set(gcf, 'Units','pixels');
   axesm('mapprojection', 'robinson', ...
         'maplatlimit', latlim, 'maplonlimit', lonlim, ...
         'grid', 'on', 'frame', 'on', 'flinewidth', 1, ...
         'parallellabel', 'off', 'meridianlabel', 'off', ...
         'MLineLocation', 60, 'PLineLocation', 20, ...
         'MLabelParallel', 'south', ...
         'labelformat', 'compass')

   surfm(glat, glon, gtmp)
   S = shaperead('landareas','UseGeoCoords',true);

   ch = colorbar('southoutside');
   p=get(ch,'position');
   set(ch,'position',p + [0 -0.04 0 -0.3*p(4)]);
   tightmap

   if do_opts
      if isfield(opts,'color');
         geoshow([S.Lat], [S.Lon],'Color',opts.color);
      end
      if isfield(opts,'title');
         th=title(opts.title);
         if isfield(opts,'titlesize')
            set(th, 'FontSize', opts.titlesize);
         else
            set(th, 'FontSize', 12)
         end
      end
      if isfield(opts,'caxis');
         caxis(opts.caxis);
      end
      if isfield(opts,'cmap');
         colormap(opts.cmap);
      end
   else
      geoshow([S.Lat], [S.Lon],'Color','k');
      colormap(jet);
   end
end


