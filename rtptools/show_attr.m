function [] = show_attr( xattr );

% function [] = show_attr( xattr );
%
% Shows the attribute text by printing each strings to the screen.
%
% Input:
%    xattr = RTP hattr or pattr (cell array with 1x3 elements per entry)
%
% Output: (none)
%

% Created 4 January 2002, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nattr=length(xattr);

disp(' ')
for ia=1:nattr
   msg1=['%%%%%%%%%%%%%%%%%%% ' char( xattr{1}(1) ) ' ' int2str(ia) ...
        ' %%%%%%%%%%%%%%%%%%%'];
   msg2=char( xattr{ia}(2) );
   msg3=char( xattr{ia}(3) );
   disp(msg1);
   disp([msg2 ' | ' msg3]);
   disp(' ')
end

%%% end of function %%%
