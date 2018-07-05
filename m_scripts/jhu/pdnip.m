function pdnip(natoms,rho,config_m);

%SH 7-29-15
%displays the x y z configuration file specifically for the Pd Ni P,
%40-40-20 ratio
    
%load config_m.txt;


halfedge = ((natoms/rho)^(1/3))/2;
%halfedge = ((4000/0.0783)^(1/3))/2;

radpd= 1.41/halfedge;
radni = 1.28/halfedge;
radp= 1/halfedge;

var = natoms/5;
con=config_m;
con=con+1;
sphereview(2, radpd, con(1:var*2,:), .9);
sphereview(2, radni, con(var*2+1:var*4,:), .8);
sphereview(2, radp, con(var*4+1:var*5,:), .6);
axis equal;
camlight right;

