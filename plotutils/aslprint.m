function [] = aslprint(fn);

% Filter for export_fig.  Seems it can't handle ~ for .pdf files.
% L. Strow, Aug. 17, 2013.
%
% Usage: llsprint(filename.ext)
%   Be sure to include either .pdf or .png extension (.ext)

% Path to export_fig
%addpath /asl/matlib/fileexchange/export_fig
addpath /asl/matlib/fileexchange/export_fig

% First get username (works on Mac too), space at end of command?
[stat username] = system('id -u -n');
username = username(1:(end-1));

% Replace ~ with appropriate full path
if length(fn) > 1
   if fn(1) == '~'
         fn = strrep(fn,'~',['/home/' username]);
   end
end

% Make background white for printing, then reset to Matlab default color
% Put the figure back into the docked/undocked state it started in
dockstate = get(gcf,'windowstyle')
% Temporarily make it a normal window so we can change it's size
set(gcf,'windowstyle','normal')
% Want figure background transparent, looks nicer in beamer slides
set(gcf,'color','none');
%  Set figure to default size
p = get(gcf,'position');
set(gcf,'position',[p(1) p(2) 560 420]);
% Export both .pdf and .png
if fn(end-3:end) == '.png'
   fn = fn(1:end-4);
elseif  fn(end-3:end) == '.pdf'
   fn = fn(1:end-4);
end
% export figure
export_fig([fn '.png'],'-m2','-transparent');
export_fig([fn '.pdf'],'-m2','-transparent');
% export_fig([fn '.png'],'-m2');
% export_fig([fn '.pdf'],'-m2');
% Set the background color back to Matlab default
set(gcf,'color',[0.8 0.8 0.8]);
% Set the dockstate back to what it was
set(gcf,'windowstyle',dockstate)
fnfig = fn;
hgsave(gcf,fnfig);
