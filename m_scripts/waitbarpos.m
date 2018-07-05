function h = waitbarpos(x, caption, left, bottom)

%SH 12-1-08
%
%calls waitbar but allows you to place it somewhere on the screen.  if left
%and bottom are unspecified, then it will be placed in the bottom left corner
%of the screen

if nargin<3
    left = 0;
    bottom = 0;
end

h = waitbar(x, caption, 'Position', [left bottom 380 50]);
