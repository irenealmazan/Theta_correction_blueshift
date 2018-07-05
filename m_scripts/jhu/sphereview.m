function h=sphereview(edge,radius,atompos,colorfrac);

% SH 6-16-05
% gives hard-sphere view of orientation of one type of atom accordign to
% atompos in a square cell defined by edge
%
% edge - size of model cube edge
% atompos - n X 3 matrix containing x y z coordinates for n atoms of same
%           type
% colorfrac - a number between 0 and 1 to which the colors scale
% radius - radius of the atom in angstroms

colormap('default');
[x,y,z]=sphere(60);
%fvc=surf2patch(x*.001,y*.001,z*.001,ones(size(z)));
%patch(fvc); 
%shading flat;

for i=1:size(atompos,1)
    xtemp=x*radius+atompos(i,1);
    ytemp=y*radius+atompos(i,2);
    ztemp=z*radius+atompos(i,3);
    h=surf(xtemp,ytemp,ztemp,'EdgeColor','none','FaceColor',colorfrac);
    %h=surf(xtemp,ytemp,ztemp,'EdgeColor','none');
    
    %fvc=surf2patch(xtemp,ytemp,ztemp,ones(size(z))*colorfrac);
    %patch(fvc);
end
    
    %shading flat;
    %axis([0 edge 0 edge 0 edge]);
    

