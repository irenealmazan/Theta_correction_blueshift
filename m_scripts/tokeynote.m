function tokeynote(fh)

if nargin<1
    fh = gcf;
end

pos = get(fh, 'position');
set(fh, 'position', [pos(1:2) 400 350]);
set(fh, 'paperpositionmode', 'auto');
set(fh,'color','w');

hline = findobj(gca, 'type','line');
w_o = get(hline,'linewidth');
set(hline, 'linewidth', 1.5);

htext = findall(findobj(gcf), 'type','text');
font_o = get(htext, 'FontName');
size_o = get(htext, 'fontsize');
set(htext, 'FontName', 'Arial');
set(htext, 'Fontsize', 14);

axsize_o = get(gca, 'fontsize');
set(gca, 'fontsize', 13);

print('-depsc', '~/temp.eps');
unix('automator ~/m_scripts/tokeynote.app');

set(gcf, 'position', pos);

for ii=1:length(hline) set(hline(ii), 'linewidth', w_o{ii}); end
for ii=1:length(htext) 
    set(htext(ii), 'fontname', font_o{ii}); 
    set(htext(ii), 'fontsize', size_o{ii});     
end
set(gca, 'fontsize', axsize_o);