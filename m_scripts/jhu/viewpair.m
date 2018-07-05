function h=viewpair(pair, atom, vertices, confile, frac, light)

nnvisual(pair.cnaidlocal, atom, vertices, confile, frac, 0, 0,1);
nnvisual(pair.atom(1).localid, atom, vertices, confile, frac, 0,0,0);
nnvisual(pair.atom(2).localid, atom, vertices, confile, frac, light, 0,0);
display(['cna indices are ' num2str(pair.cna)]);
