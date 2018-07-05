function [pairs bh] = cnalocal(atom,vertices,atden,nelem,radii,r)

natoms = length(atom);
halfedge = ((1/atden)*natoms)^(1/3) /2
radii = radii/halfedge;
paircnt = 1;
rspacing = r(2)-r(1);
bh=zeros(length(r),1);
bh2=zeros(length(r),1);
pairidmatrix=[];

for i=1:natoms %go through all the atoms as an origin
    
    if mod(i,5)==0 display(['calculating origin atom ' num2str(i)]);
    end

    %now determine the nearest neighbor list
    clear localnn originnn;
    clear templocalnn;
    
    localnn = atom(i).neighborlist;
    originnn= atom(i).neighborlist;
    
    for j=1:atom(i).nn
        tempindex = originnn(j);
        templocalnn = atom(tempindex).neighborlist;
        localnn = horzcat(localnn, templocalnn);
    end
   
    localnn = unique(localnn)
    nlocalnn = length(localnn);
    
    for j=1:nlocalnn %go through all the columns in this row without being redundant
        
        tempindex = localnn(j);
        
        pairid = sqrt(atom(i).pos(1)) * sqrt(atom(tempindex).pos(1));
        
        if ~ismember(pairid, pairidmatrix)
            
            delx = atom(i).pos(1) - atom(tempindex).pos(1);
            dely = atom(i).pos(2) - atom(tempindex).pos(2);
            delz = atom(i).pos(3) - atom(tempindex).pos(3);

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

            if i~=tempindex  
                
                nni = atom(i).nn;
                nnj = atom(tempindex).nn;
                
                tempindex
                pause;
                
                nnidi = atom(i).neighborlist;
                nnidj = atom(tempindex).neighborlist;

                clear cnaid;
                cnaid = intersect(nnidi,nnidj); %gives a list of intersecting in ascending order

                %display(num2str([i tempindex dist]));


                if ~isempty(cnaid) %then this pair is a candidate for cna analysis

                    %[i j dist*halfedge]
                    %cnaid

                    %fill in the second and third cna indices
                    if length(cnaid) >1

                        if bin<=length(r) & dist>0
                            bh2(bin) = bh2(bin)+1;
                        end

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

                        [L,tree]= lindex(atom, cnaid);

                        cnaresult(1) = length(cnaid); %this is the first cna index
                        %it records the number of 
                        cnaresult(2) = length(unique(bonds));
                        cnaresult(3) = max(L);
                    end  %end if length(cnaid)>1         


                    if length(cnaid) ==1
                        cnaresult = [1 0 0];
                    end

                    pairs(paircnt).atom(1) = atom(i);
                    pairs(paircnt).atom(2) = atom(tempindex);
                    pairs(paircnt).cnaid = cnaid;
                    pairs(paircnt).cna = cnaresult;
                    pairs(paircnt).dist = dist;

                    pairs(paircnt).centroid(1) = (atom(i).pos(1)+ ...
                        atom(tempindex).pos(1)) / 2;
                    pairs(paircnt).centroid(2) = (atom(i).pos(2)+ ...
                        atom(tempindex).pos(2)) / 2;
                    pairs(paircnt).centroid(3) = (atom(i).pos(3)+ ...
                        atom(tempindex).pos(3)) / 2;
                    
                    pairidmatrix(paircnt)=pairid;
                    
                    paircnt = paircnt + 1;

                    %[i j dist*halfedge cnaresult];
                    %cnaid;

                    %pause(.05)
                end %end if ~isempty(cnaid) 
        
            end %if i ~=j
        
        end %if ~ismember(pairid, pairidmatrix)
    end %end for j
   
end %end for i
bh(:,2)=bh2;
plot(r,bh)
    
    