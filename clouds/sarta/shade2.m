function haha = shade2(fig,X0,Y0,W,H,color,transperancy)

%% function haha = shade(fig,x0,y0,w,h,color,transperancy)
%% shades figure between [X0 X0+W X0+W X0 X0],[Y0 Y0 Y0+H Y0+H Y0]
%% using "color" and transperancy [0 0.5 1] = transperant inbetween opaque 
%%
%% XO = xstart, W = width
%% Y0 = ystart, H = height

figure(fig)
hh = patch([X0 X0+W X0+W X0 X0],[Y0 Y0 Y0+H Y0+H Y0],[0 0 0 0 0],... 
           'FaceColor',color,'Edgecolor',color); 
set(hh,'FaceAlpha',transperancy);  %%[0 0.5 1 = transperant inbetween opaque 
set(hh,'EdgeAlpha',transperancy);  %%[0 0.5 1 = transperant inbetween opaque 
