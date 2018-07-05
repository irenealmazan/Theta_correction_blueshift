function [posincol num] = columnpos(pos, x, y, Q, edge)

halfedge= edge/2;
radius = .61/Q;
totat= length(pos);
posincol=[];

display(['halfedge is ' num2str(halfedge) ', radius is ' num2str(radius)]);

plot3(pos(:,1),pos(:,2), pos(:,3),'o', x,y,0, 'r*');

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

hold; plot3(posincol(:,1), posincol(:,2), posincol(:,3),'g*'); hold off

num=length(posincol);

