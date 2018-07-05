clear v
clear vstat
clear m
clear vl

%bort = 'bmp';
bort = 'tif';
suffix = '';
num = 5;

A=double(imread(['A' suffix '.' bort]));
B=double(imread(['B' suffix '.' bort]));
C=double(imread(['C' suffix '.' bort]));
D=double(imread(['D' suffix '.' bort]));
if num>4
    E=double(imread(['E' suffix '.' bort]));
end
if num>5
    F=double(imread(['F' suffix '.' bort]));
end

%  A=outlyers(A,4);
%  B=outlyers(B,4);
%  C=outlyers(C,4);
%  D=outlyers(D,4);
%  E=outlyers(E,4);
%  if num ==6
%      
%      F=outlyers(F,4);
%  end


[v(1) m(1) vstat(1) vl(1)]= expvar(A);
[v(2) m(2) vstat(2) vl(2)]= expvar(B);
[v(3) m(3) vstat(3) vl(3)]= expvar(C);
[v(4) m(4) vstat(4) vl(4)]= expvar(D);
if num > 4
    [v(5) m(5) vstat(5) vl(5)]= expvar(E);
end
if num >5
    [v(6) m(6) vstat(6) vl(6)]= expvar(F);
end

g=[.2261 .4522 .6219 .7915 .9611 1.1307];

plot(g(1:length(v)), v, g(1:length(v)), vstat,'o', g(1:length(v)), vl,'+');