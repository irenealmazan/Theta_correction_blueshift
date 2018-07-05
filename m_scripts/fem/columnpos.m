function [posincol num] = columnpos(pos, x, y, colrad, edge)
% Calculate the positiongs of atoms in a column projected 
% onto the x,y plane, relative to the center of the column
% Input:
% pos = array of atomic positions in model [x y z zed] (zed=atomic number)
% x,y = center of column
% colrad = radius of column (=0.61/Q where Q=obj aperture radius)
% edge = box size (Angstroms)
%
% Output:
% posincol = array of positions of atoms in desired column
% num = number of atoms in column

% Stephan Hruszkewycz and Todd Hufnagel (hufnagel@jhu.edu)
% 28-Mar-2009

halfedge= edge/2;
radius = colrad;
totat= length(pos);
posincol=[];

%display(['halfedge is ' num2str(halfedge) ', radius is ' num2str(radius)]);

%plot3(pos(:,1),pos(:,2), pos(:,3),'o', x,y,0, 'r*');

% Calculate the projected positions, taking account of
% periodic boundary conditions
for i=1:totat
    
    delx = pos(i,1)-x;
    dely = pos(i,2)-y;
    
    if abs(delx) > halfedge
        delx = edge - abs(delx);
    end
    if abs(dely) > halfedge
        dely = edge - abs(dely);
    end
    
    sigma = sqrt( delx.^2 + dely.^2 );
    if sigma <= radius
        posincol = [posincol; i];
    end
    
end
posincol = pos(posincol, [1 2 4]);
num=length(posincol);

hold off;
%plot3(posincol(:,1), posincol(:,2), posincol(:,3),'go');
plot(posincol(:,1),posincol(:,2),'o');
axis([-edge/2 edge/2 -edge/2 edge/2]);
pause(.1)
title(['column position: ' num2str(x) num2str(y)]);

% Make the projected positions relative to the center of the column,
% taking into account the periodic boundary conditions
for i=1:num
	delx = posincol(i,1) -x;
	dely = posincol(i,2) -y;

    if abs(delx) <= halfedge
        posincol(i,1) = posincol(i,1) - x;
    end
    if abs(delx) > halfedge
        posincol(i,1) = posincol(i,1) - sign(delx)*edge - x;
    end

    if abs(dely) <= halfedge
        posincol(i,2) = posincol(i,2) - y;
    end
    if abs(dely) > halfedge
        posincol(i,2) = posincol(i,2) - sign(dely)*edge - y;
    end

end
hold on; plot(posincol(:,1),posincol(:,2),'ro');
axis([-edge/2 edge/2 -edge/2 edge/2]);
pause(.1); 
hold off

