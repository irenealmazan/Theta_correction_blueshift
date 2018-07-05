function sup = shrinkwrap(obj, sigma, isofrac)

%SH 11-11-08
%returns a new shrink-wrapped support defined by convovling the object
%(obj) with a 3d normalized gaussian function with half-width defined by
%sigma.  an isosurface of the blurred object is extracted at an isovalue of
% (isofrac)*(maximum value in obj).  to ensure that internal holes or
% divots in the obj are not taken as zero, only the convex hull of the
% isosurface points is used to make the new support (sup).

obj=convolvegauss(obj.object, sigma);

%using vertices from isosurface
x=[1:size(obj,2)];
y=[1:size(obj,1)];
z=[1:size(obj,3)];

[X,Y,Z]= meshgrid(x,y,z);

scale= max(max(max(abs(obj))));
scale= scale * isofrac;

[f,v]= isosurface(X,Y,Z,abs(obj),scale);

K = convhulln(v);
vcon= v(unique(K), :);
vcon = [ceil(vcon) ; floor(vcon)];
vcon = unique(vcon, 'rows');


% clf
% p=patch('Faces',f,'Vertices',v);
% set(p,'FaceColor','red','EdgeColor','black');
% camlight 
% alpha(.3)
% axis equal
% hold
% plot3(vcon(:,1),vcon(:,2), vcon(:,3), 'c.'); pause



%add points to vcon so that the delaunay calculation is more accurate and
%does not need to be done twice

minz = floor(min(vcon(:,3)));
maxz = ceil(max(vcon(:,3)));
minrow = floor(min(vcon(:,1)));
maxrow = ceil(max(vcon(:,1)));
mincol = floor(min(vcon(:,2)));
maxcol = ceil(max(vcon(:,2)));

% counter =1;
% addvcon =[];
% 
% for zaxis = minz:maxz
%     zaxis
%     inplane = vcon(find(vcon(:,3)==zaxis), :);
%     
%     for i= min(inplane(:,2)): max(inplane(:,2)) %cols
% 
%         coltemp = [min(inplane(:,1)): max(inplane(:,1))];
%         coltemp = [coltemp' (ones(size(coltemp))*i)' (ones(size(coltemp))*zaxis)']
%         setdiff(coltemp, inplane(find(inplane(:,2)==i),:), 'rows')
%         
%         addvcon = [addvcon; setdiff(coltemp, inplane, 'rows')];
%         clear coltemp
%     end
% 
%     plot(inplane(:,1), inplane(:,2),'o'); hold on;
%     
%     if addvcon & length(find(addvcon(:,3)==zaxis))
%         plot(addvcon(find(addvcon(:,3)==zaxis),1), addvcon(find(addvcon(:,3)==zaxis),2),'ro');
%     end
% 
%     hold off
%     pause(.2)
% 
% end
% 
% pause


%use the points in vcon to 

T=delaunayn(vcon);
counter =1;
clear pts

for zaxis=minz: maxz
    
    for col = mincol:maxcol
        
       for row = minrow:maxrow
           
           pts(counter,:) = [row, col, zaxis];
           counter=counter+1;
           
       end
        
    end
        
end


% p=tetramesh(T,vcon);
% set(p,'FaceColor','red','EdgeColor','none');
% camlight right; 
% pause


k=dsearchn(vcon, T, pts, -1);
pts = pts(find(k ~= -1),:);


%do it again to pick up missed internal points

T=delaunayn(pts);
counter =1;
clear pts2

for zaxis=minz: maxz
    
    for col = mincol:maxcol
        
       for row = minrow:maxrow
           
           pts2(counter,:) = [row, col, zaxis];
           counter=counter+1;
           
       end
        
    end
        
end

k=dsearchn(pts, T, pts2, -1);
pts = pts2(find(k ~= -1),:);


%make the support
sup = zeros(size(obj));

for i=1:length(pts)
    sup(pts(i,2), pts(i,1), pts(i,3))=1;
end
