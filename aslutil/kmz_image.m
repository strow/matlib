function kmz_image(filename)
% function kmz_image(filename)
%
% Create a kmz file with the image created by a mapping routine such as simplemap.
%
% INPUTS:
%  gca (MATLAB Figure)  - raw graphic image to export to a kmz file
%  filename             - output file name (must end in .kmz)
%
% OUTPUT:
%  a file is written to the file system
%
% EXAMPLE:
%  simplemap(rand(10)*10,rand(10)*10,rand(10))
%  %simplemap(p.rlat,p.rlon,p.robs1(1291,:),0.15,'auto','ngrow') % airs data
%  kmz_image('test.kmz')  % write out the image

% Written by Paul Schou, March 2013

data = [];
for g = get(gca,'Children')'
  if isfield(get(g),'CData')
    data = get(g,'CData');
    alpha = get(g,'AlphaData');
    yd = get(g,'YData'); yd = yd([1 end]);
    yd = yd + [-1 1] * diff(yd)  / (size(data,2)-1) / 2;
    xd = get(g,'XData'); xd = xd([1 end]);
    xd = xd + [-1 1] * diff(xd)  / (size(data,1)-1) / 2;
    break
  end
end
if isempty(data)
  error('No image found on current figure handle')
end

if ~strcmp(filename(max(1,end-3):end),'.kmz')
  error('Filename must end in .kmz')
end

bn = filename(1:end-4);
kml_fn = [bn '.kml'];
png_fn = [bn '.png'];

cm = colormap;
cbins = linspace(min(caxis),max(caxis),size(cm,1));
[x bins] = histc(data,cbins);
bins(data >= max(cbins)) = length(cbins)-1;
bins(data <= min(cbins)) = 1;
bins(bins < 1) = 1;
%imagesc(reshape(cm(bins,:),[size(data) 3]));
%keyboard
%alpha = fliplr(alpha);
if strcmp(get(gca,'YDir'),'normal')
  bins = flipud(bins);
  alpha = flipud(alpha);
end
imwrite(reshape(cm((bins),:),[size(data) 3]),png_fn,'Alpha',alpha+0); 

% One can fine documentation about the xml + kml details here:
%  https://developers.google.com/kml/documentation/kml_tut
fid = fopen(kml_fn,'w');
fwrite(fid, [...
  '<?xml version="1.0" encoding="UTF-8"?>' 10 ...
  '<kml xmlns="http://www.opengis.net/kml/2.2">' 10 ...
  '<GroundOverlay>' 10 ...
  '	<name>data</name>' 10 ...
  ' <Icon>' 10 ...
  '   <href>' png_fn '</href>' 10 ...
  '  </Icon>' 10 ...
  '  <LatLonBox>' 10 ...
  '    <north>' num2str(max(yd)) '</north>' 10 ...
  '    <south>' num2str(min(yd)) '</south>' 10 ...
  '    <east>' num2str(max(xd)) '</east>' 10 ...
  '    <west>' num2str(min(xd)) '</west>' 10 ...
  '    <rotation>0</rotation>' 10 ...
  '  </LatLonBox>' 10 ...
  '  <LookAt>' 10 ...
  '    <longitude>' num2str(mean(xlim)) '</longitude>' 10 ...
  '    <latitude>' num2str(mean(ylim)) '</latitude>' 10 ...
  '    <altitude>' num2str(atan(max(diff(ylim),diff(xlim))/30)/pi*20000000) '</altitude>' 10 ...
  '    <altitudeMode>clampToGround</altitudeMode>' 10 ...
  '  </LookAt>' 10 ...
  '</GroundOverlay>' 10 ...
  '</kml>']);
  %'    <altitude>' num2str(atan(max(10)/35)/pi*20000000) '</altitude>' 10 ...

fclose(fid);

zip([bn '.zip'],{kml_fn png_fn})
movefile([bn '.zip'], filename)
delete(kml_fn)
delete(png_fn)
