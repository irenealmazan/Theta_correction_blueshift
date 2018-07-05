function [pairs bh] = cna(atom,vertices,atden,nelem,radii,r,natoms,natomstot)

halfedge = ((1/atden)*natomstot)^(1/3) /2
radii = radii/halfedge;
paircnt = 1;
rspacing = r(2)-r(1);
bh=zeros(length(r),1);
h = waitbar(0,'calculating pairs');


for i=1:natoms %go through all the rows
    
    %if mod(i,25)==0 display(['calculating row ' num2str(i)]);
    %end
%     %first do the voronoi index
%     atom(i).vorind = [0 0 0 0 0 0 0 0 0]; % n3 n4 n5 n6 n7 n8, or number of nn that create n sided polygons
%         %this is a new field added to atom.  edges are screwed up the small
%         %faces left on the polygons.  have to re-calculate the voronoi here
%         %with only the right amount of nn.
%         
%     nn = atom(i).nn;
%     
%     for j=1:nn
%         poledge = atom(i).neighbor(j).nedges;  %store how many edges this face has
%         atom(i).vorind(poledge-2) = atom(i).vorind(poledge-2)+1;  %increment the appropriate bin
%     end
    
    for j=i:natoms %go through all the columns in this row without being redundant
        
        delx = atom(i).pos(1) - atom(j).pos(1);
        dely = atom(i).pos(2) - atom(j).pos(2);
        delz = atom(i).pos(3) - atom(j).pos(3);
            
        %check bc in x direction
        if delx >= 1
            delx = 2-delx;
        end
        if delx < -1
            delx = 2+delx;
        end

        %check bc in y direction
        if dely >= 1
            dely = 2-dely;
        end
        if dely < -1
            dely = 2+dely;
        end

        %check bc in z direction
        if delz >= 1
            delz= 2-delz;
        end
        if delz < -1
            delz = 2+delz;
        end
        
        dist = sqrt(delx^2 + dely^2 + delz^2);
        
        bin = round( dist*halfedge / rspacing) +1;
    
        if bin<=length(r) & dist>0
            bh(bin) = bh(bin)+1;
        end
        
        %if  bin==48 display(num2str([i j])); end
        
        if dist < 5*max(radii) & i~=j  %then they are close enough to do cna
        
            nni = atom(i).nn;
            nnj = atom(j).nn;
            
            clear nnidi nnidj
            nnidi = zeros(nni,1);
            nnidj = zeros(nnj,1);
            nnlocalidi = zeros(nni,1);
            nnlocalidj = zeros(nnj,1);
            
            for k=1:nni
                nnidi(k) = atom(i).neighbor(k).id;
                nnlocalidi(k) = atom(i).neighbor(k).localid;
            end
            for k=1:nnj
                nnidj(k) = atom(j).neighbor(k).id;
                nnlocalidj(k) =atom(j).neighbor(k).localid;
            end
            
            clear cnaid cnaidlocal;
            cnaid = intersect(nnidi,nnidj); %gives a list of intersecting in ascending order
            cnaidlocal= intersect(nnlocalidi, nnlocalidj);
            
            if ~isempty(cnaid) %then this pair is a candidate for cna analysis
                
                %[i j dist*halfedge]
                %cnaid
                
                %fill in the second and third cna indices
                if length(cnaid) >1

                    clear nnlisttemp intersecttemp bonds;
                    nnlisttemp = zeros(length(cnaid),30);
                    intersecttemp = zeros(length(cnaid), 10);
                    bonds = zeros(200,1); %tons of room just in case

                    kcounter=0;
                    for l=1:length(cnaid)

                        atindexl = cnaid(l);
                        nnl = atom(atindexl).nn;

                        for m=1:nnl %fill in the nnlisttemp matrix with the neigbors of all cna atoms
                            nnlisttemp(l,m) = atom(atindexl).neighbor(m).id;
                        end

                        intersectl = intersect(nnlisttemp(l,:), cnaid);

                        for m =1:length(intersectl)
                            kcounter=kcounter+1;
                            bonds(kcounter)= sqrt(intersectl(m)) * sqrt(atindexl);
                                %this gives a fairly unique number that
                                %describes this pair of atomic indices. 
                        end


                    end    

                    bonds = bonds(find(bonds~=0));
                    
                    %[L,tree]= lindex(atom, cnaid);
                    [L,tree]= lindex(atom, cnaidlocal);
                    
                    cnaresult(1) = length(cnaid); %this is the first cna index
                    %it records the number of 
                    cnaresult(2) = length(unique(bonds));
                    cnaresult(3) = max(L);
                end  %end if length(cnaid)>1         
                
                
                if length(cnaid) ==1
                    cnaresult = [1 0 0];
                end
                
                pairs(paircnt).atom(1) = atom(i);
                pairs(paircnt).atom(2) = atom(j);
                pairs(paircnt).cnaid = cnaid;
                pairs(paircnt).cnaidlocal = cnaidlocal;
                pairs(paircnt).cna = cnaresult;
                pairs(paircnt).dist = dist;

                pairs(paircnt).centroid(1) = (atom(i).pos(1)+ ...
                    atom(j).pos(1)) / 2;
                pairs(paircnt).centroid(2) = (atom(i).pos(2)+ ...
                    atom(j).pos(2)) / 2;
                pairs(paircnt).centroid(3) = (atom(i).pos(3)+ ...
                    atom(j).pos(3)) / 2;
                
                
                paircnt = paircnt + 1;
                
                %[i j dist*halfedge cnaresult];
                %cnaid;
                
                %pause(.05)
            end %end if ~isempty(cnaid) 
            
            
        end %end if dist <
        
        
    end %end for j
    waitbar(i/natoms);
end %end for i

close(h);
plot(r,bh)
    
    