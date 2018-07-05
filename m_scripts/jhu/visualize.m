function visualize(config);

%SH 7-29-15
%displays the x y z configuration file, giving different colors and
%proportional sizes to the different atoms according to the configuration
%file given
%
%input :  config - string of the file name containing atomic coordinates
%that were output from the rmc simulation 
    

%first read the config file
fid1=fopen(config, 'r');
rho = str2num(fgetl(fid1));
atomtypes = str2num(fgetl(fid1));
numdifat=atomtypes;
color=[.9 .5 0; .9 .9 0; 0 0 .5 ; .3 .8 .8 ; .5 .1 .1 ; .3 .9 .3];

for i=1:atomtypes
    row = str2num(fgetl(fid1));
    numdifat(i) = row(1);
    radius(i)= row(2);
end

con = fscanf(fid1, '%f %f %f', [3 inf]);
con = con';

%con=con+1; %this brings the atoms to a 0 to 2 scale

natoms=sum(numdifat);

%put into box coordinates
halfedge = ((natoms/rho)^(1/3))/2;

radius=radius/halfedge;
counter=1;
hold on;

for i=1:atomtypes %goes through each atom type
    for j=1:numdifat(i) %now go through each atom of this type
        h=sphereview(2,radius(i),con(counter,:),color(i,:));
        %if i ~= 3
            %display('hi')
        %    set(h,'FaceAlpha',.3);  %this makes this atom transparent
        %end  
        counter=counter+1;
    end
end

hold off;
axis equal;
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);
camlight right;
lighting phong;
