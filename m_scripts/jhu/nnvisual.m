function nnvisual(which2draw, atom, vertices, confilename, atsizefrac, light, shownn, showvoron)


%first read the config file
fid1=fopen(confilename, 'r');
rho = str2num(fgetl(fid1));
atomtypes = str2num(fgetl(fid1));
numdifat=atomtypes;
color=[.9 .5 0; .9 .9 0; 0 0 .5 ; .3 .8 .8 ; .5 .1 .1 ; .3 .9 .3];

for i=1:atomtypes
    row = str2num(fgetl(fid1));
    numdifat(i) = row(1);
    radius(i)= row(2);
end

config = fscanf(fid1, '%f %f %f', [3 inf]);
config = config';

fclose(fid1);

% for i=1:sum(numdifat)
%     config(i, :) = str2num(fgetl(fid1));
% end


%begin by drawing the central atom and its nearest neighbors of each 
for atomcnt=1:length(which2draw)
    whichatom = which2draw(atomcnt);

    hold on;
    halfedge = ((sum(numdifat)/rho)^(1/3))/2;
    type = atom(whichatom).type;
    color(type,:);
    h=sphereview(2, radius(type)*atsizefrac/halfedge, atom(whichatom).pos, color(type,:));

    if(shownn)
        for i=1:atom(whichatom).nn 

            x = atom(whichatom).pos(1) + atom(whichatom).neighbor(i).dpos(1);
            y = atom(whichatom).pos(2) + atom(whichatom).neighbor(i).dpos(2);
            z = atom(whichatom).pos(3) + atom(whichatom).neighbor(i).dpos(3);

            type = atom(whichatom).neighbor(i).type;

            h=sphereview(2, radius(type)*atsizefrac/(halfedge), [x y z], color(type,:));

        end
    end

    vertind = atom(whichatom).vertices;
    vertind = vertind(find(vertind>0));
    vertexlist = vertices(vertind,:);
    vertexlist = vertexlist;

    %draw the polyhedron

    if showvoron
        X=[vertexlist(:,1) vertexlist(:,2) vertexlist(:,3)];
        K=convhulln(X);

        dhull=[1 2 3 1]; %picks out the trinangle vertices and completes the triangle

        for k=1:size(K,1)
            l=K(k,dhull);
            h(i)=patch(X(l,1),X(l,2),X(l,3),k,'EdgeColor','none','FaceColor','g','FaceAlpha',.6);
        end
        %set(findobj('Type','patch'),'EdgeColor','none')
        %set(findobj('Type','patch'),'FaceColor','g')
    end
end
    
axis equal;

if(light)
    camlight right; 
    lighting phong;
end

hold off;   



    


       