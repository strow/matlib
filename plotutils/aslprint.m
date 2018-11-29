function [] = aslprint(varargin);

% Filter for export_fig.
%
% Usage: aslprint(filename,anyvar)
%   Do not include filename extension!
%   If anyvar exists, only png output
%  
nargs = length(varargin);
fn = varargin{1}
if nargs == 2
   save_png_only = true;
else
   save_png_only = false;
end

% Path to export_fig
addpath /asl/matlib/fileexchange/export_fig

% First get username (works on Mac too), space at end of command?
[stat username] = system('id -u -n');
username = username(1:(end-1));

% Replace ~ with appropriate full path
if length(fn) > 1
   if fn(1) == '~'
      switch computer
         case 'MACI64'
         fn = strrep(fn,'~',['/Users/' username]);
         case 'GLNXA64'
         fn = strrep(fn,'~',['/home/' username]);
        otherwise
          disp('Unknown matlab computer type?');
      end
   end
end

% Make background white for printing, then reset to Matlab default color
% Put the figure back into the docked/undocked state it started in
dockstate = get(gcf,'windowstyle');
% Temporarily make it a normal window so we can change it's size
set(gcf,'windowstyle','normal')
% Want figure background transparent, looks nicer in beamer slides
set(gcf,'color','none');
%set(gcf,'color',[1 1 1]);
%  Set figure to default size
p = get(gcf,'position');
set(gcf,'position',[p(1) p(2) 560 420]);
% Export both .pdf and .png
if fn(end-3:end) == '.png'
   fn = fn(1:end-4);
elseif  fn(end-3:end) == '.pdf'
   fn = fn(1:end-4);
end

warning('off');

% export figure
export_fig([fn '.png'],'-m2','-transparent');
if ~save_png_only
   export_fig([fn '.pdf'],'-m2','-transparent');
end

% Set the background color back to Matlab default
set(gcf,'color',[0.94 0.94 0.94]);
% Set the dockstate back to what it was
set(gcf,'windowstyle',dockstate)
fnfig = fn;
hgsave(gcf,fnfig);

warning('on');

