function llsprint(varargin);
% Usage
% L. Strow, Jan. 17, 2015
%

nargs = length(varargin);
fn = varargin{1};
if nargs == 1
   aslprint(fn);
elseif nargs == 2
   ovar = varargin{2};
   aslprint(fn,ovar);
end

fnpng = [fn '.png'];

eval_str = ['scp -p ' fnpng ' spud:Dropbox'];

unix(eval_str);

