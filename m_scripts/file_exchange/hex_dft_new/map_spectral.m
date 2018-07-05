function out = map_spectral(hex_add)
global N;
global C;
digits = fliplr(arrayfun(@(x) str2double(x),num2str(hex_add)));
Tf = (1/3) * [-1 2 -1;2 -1 -1];
h = arrayfun(@(x) N{x}.' *Tf*C{digits(x)+1},1:length(digits),'UniformOutput',0);
if length(digits)==1
    out = h{1};
else
    out = sum(cat(3,h{:}),3);
end
