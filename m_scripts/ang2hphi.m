function [HKL]=ang2hphi(nu, chi, eta, del, dspacing)

cn = cosd(nu);
sn = sind(nu);

cd = cosd(del);
sd = sind(del);

cc = cosd(chi);
sc = sind(chi);

ce = cosd(eta);
se = sind(eta);

H = (ce - cd * cn * ce - se * sd) * sc + cd * sn * cc;
K = -(ce - cd * cn * ce - se * sd) * cc + cd * sn * sc;
L = ce * sd + se * (1 - cn * cd); 

HKL=[H K L];

a1=[dspacing 0 0];
a2=[0 dspacing 0];
a3=[0 0 dspacing];

Rxeta=[1   0   0; 
    0   cosd(eta)  sind(eta);
    0   -sind(eta) cosd(eta)];

chi = -chi;
Rychi=[cosd(chi) 0 -sind(chi);
    0   1   0;
    sind(chi) 0 cosd(chi)];

Rxdel=[1   0   0; 
    0   cosd(del)  sind(del);
    0   -sind(del) cosd(del)];

Rynu=[cosd(nu) 0 -sind(nu);
    0   1   0;
    sind(nu) 0 cosd(nu)];

a1 = (Rychi*a1')';
a2 = (Rychi*a2')';
a3 = (Rychi*a3')';

display(['a1: ' num2str(a1)]);
display(['a2: ' num2str(a2)]);
display(['a3: ' num2str(a3)]);

a1 = (Rxeta*a1')';
a2 = (Rxeta*a2')';
a3 = (Rxeta*a3')';

display('Rxeta');
display(num2str(Rxeta)); display(' ');
display('Rychi');
display(num2str(Rychi)); display(' ');

display(['a1: ' num2str(a1)]);
display(['a2: ' num2str(a2)]);
display(['a3: ' num2str(a3)]);

det=[541852  -28212.3  840000];
beamdir=[0 0 1];
lambda = .0001239;

qvec(1) = (det(1)/norm(det) - beamdir(1))/lambda;
qvec(2) = (det(2)/norm(det) - beamdir(2))/lambda;
qvec(3) = (det(3)/norm(det) - beamdir(3))/lambda;

display(['qvec: ' num2str(qvec)]);

H=dot(qvec, a1);
K=dot(qvec, a3);
L=dot(qvec, a2);

HKL=[H K L];


display(' ');
vec=[.4 10 3];
mat=[3 5 6;
    9 3 30;
    11 80 8];
display(num2str( (mat*vec')'));

vecans(1) = dot(vec, mat(1,:));
vecans(2) = dot(vec, mat(2,:));
vecans(3) = dot(vec, mat(3,:));
display(num2str(vecans));

