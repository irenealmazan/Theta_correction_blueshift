function [qtr rtr] = pixelconvert(params, sizearray)

%SH 11-14-08
%
%params(1) = 

gamma = params.gamma;
delta = params.delta;
lambda = params.lambda;
arm = params.arm;
dpx = params.dpx;
dpy = params.dpy;
dth = (2*pi)/360 * params.dth;

dqdpx(1) = -cosd(delta) * cosd(gamma);
dqdpx(2) = 0;
dqdpx(3) = sind(delta) * cosd(gamma);

dqdpy(1) = sind(delta) * sind(gamma);
dqdpy(2) = -cosd(gamma);
dqdpy(3) = cosd(delta) * sind(delta);

dqdth(1) = -cosd(delta) * cosd(gamma) +1;
dqdth(2) = 0;
dqdth(3) = sind(delta) * cosd(gamma);

Astar(1) = (2*pi/lambda)*dpx*dqdpx(1);
Astar(2) = (2*pi/lambda)*dpx*dqdpx(2);
Astar(3) = (2*pi/lambda)*dpx*dqdpx(3);

Bstar(1) = (2*pi/lambda)*dpy*dqdpy(1);
Bstar(2) = (2*pi/lambda)*dpy*dqdpy(2);
Bstar(3) = (2*pi/lambda)*dpy*dqdpy(3);

Cstar(1) = (2*pi/lambda)*dth*dqdth(1);
Cstar(2) = (2*pi/lambda)*dth*dqdth(2);
Cstar(3) = (2*pi/lambda)*dth*dqdth(3);

display(['Astar is:   ' num2str(Astar)]);
display(['Bstar is:   ' num2str(Bstar)]);
display(['Cstar is:   ' num2str(Cstar)]);

denomcp = cross(Bstar, Cstar);
denom = dot(Astar, denomcp);

Axdenom = cross(Bstar, Cstar);
Bxdenom = cross(Cstar, Astar);
Cxdenom = cross(Astar, Bstar);

A(1) = 2*pi*Axdenom(1) / denom;
A(2) = 2*pi*Axdenom(2) / denom;
A(3) = 2*pi*Axdenom(3) / denom;

B(1) = 2*pi*Bxdenom(1) / denom;
B(2) = 2*pi*Bxdenom(2) / denom;
B(3) = 2*pi*Bxdenom(3) / denom;

C(1) = 2*pi*Cxdenom(1) / denom;
C(2) = 2*pi*Cxdenom(2) / denom;
C(3) = 2*pi*Cxdenom(3) / denom;

Astarmag = sqrt(dot(Astar,Astar));
Bstarmag = sqrt(dot(Bstar,Bstar));
Cstarmag = sqrt(dot(Cstar,Cstar));
Amag = sqrt(dot(A,A));
Bmag = sqrt(dot(B,B));
Cmag = sqrt(dot(C,C));

fprintf('\n');
display(['Astar magnitude ' num2str(Astarmag) ', A magnitude: ' num2str(Amag)]);
display(['Bstar magnitude ' num2str(Bstarmag) ', B magnitude: ' num2str(Bmag)]);
display(['Cstar magnitude ' num2str(Cstarmag) ', C magnitude: ' num2str(Cmag)]);

fprintf('\n');
display(['a   ' num2str(A)]);
display(['b   ' num2str(B)]);
display(['c   ' num2str(C)]);

rT(:,1) = A;
rT(:,2) = B;
rT(:,3) = C;

qT(:,1) = Astar;
qT(:,2) = Bstar;
qT(:,3) = Cstar;

[oldX, oldY, oldZ] = meshgrid(1:sizearray(2), 1:sizearray(1), 1:sizearray(3));

rtr.X = zeros(size(oldX));
rtr.Y = zeros(size(oldY));
rtr.Z = zeros(size(oldZ));
qtr.X = zeros(size(oldX));
qtr.Y = zeros(size(oldY));
qtr.Z = zeros(size(oldZ));

for i=1:sizearray(3) %z
    for j=1:sizearray(2) %cols
        for k=1:sizearray(3) %rows
            
            qtr.X(k,j,i) = sum([oldX(k,j,i) oldY(k,j,i) oldZ(k,j,i)] .* qT(1,:));
            qtr.Y(k,j,i) = sum([oldX(k,j,i) oldY(k,j,i) oldZ(k,j,i)] .* qT(2,:));
            qtr.Z(k,j,i) = sum([oldX(k,j,i) oldY(k,j,i) oldZ(k,j,i)] .* qT(3,:));

            rtr.X(k,j,i) = qtr.X(k,j,i) / sizearray(1);
            rtr.Y(k,j,i) = qtr.Y(k,j,i) / sizearray(2);
            rtr.Z(k,j,i) = qtr.Z(k,j,i) / sizearray(3);
        end
    end
end