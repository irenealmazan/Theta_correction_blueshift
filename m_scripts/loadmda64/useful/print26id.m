function print26id(scale)

if nargin<1 scale=1; end

set(gcf, 'PaperPositionMode', 'auto');
pos=get(gcf, 'PaperPosition');
set(gcf, 'PaperPosition', [.5 .5 scale*pos(3)+.5 scale*pos(4)+.5]);
print;