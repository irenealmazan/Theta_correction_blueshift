function out = map_spatial(hex_add)
global N;
global C;
digits = fliplr(arrayfun(@(x) str2double(x),num2str(hex_add))); % Store the digits in the Hexagonal addressing
Ts = (1/3) * [1 1 -2;-1 2 -1];
h = arrayfun(@(x) N{x}*Ts*C{digits(x)+1},1:length(digits),'UniformOutput',0); % Store the responses for individual digits
% Now store them as a whole (Eq 3.52)
if length(digits)==1
    out = h{1};
else
    out = sum(cat(3,h{:}),3);
end
