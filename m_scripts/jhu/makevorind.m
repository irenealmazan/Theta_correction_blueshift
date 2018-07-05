function [atom vorhist]= makevorind(atom, numdifat);

%SH commented 9-12-08
%
% this function takes an atom imput and a number of atoms as input.  the
% atom structure is returned as output with a new field 'vorind,'
% consisting of seven integers.  the first integer lists the number of 3
% edged faces on the voronoi polygon of a particular atom.  The next
% integer is the number of 4-edged faces, and so on.  
% 
% The output structure 'vorhist' is a compilation of vornonoi index
% statistics broken down by atomic species.  vorhist(1).vorhist is a
% histogram of the number of instances of a given voronoi index sorted from
% most numerous to least for species 1.  One such structure exists for each
% of n atom types e.g. vorhist(n).vorhist.  vorhist(n).id(m) is the
% voronoi index that that corresponds to the mth histogram value for atom
% n.

for i=1:length(atom)
    atom(i).vorind = [0 0 0 0 0 0 0];
    
    for j=1:atom(i).nn
        
        edgenum = atom(i).neighbor(j).nedges;
        
        if edgenum >2 & edgenum <10
            atom(i).vorind(edgenum-2) = atom(i).vorind(edgenum-2) +1;
        end
    end
end

indexcount=1;
for i=1:length(atom)
    
    if mod(i,1000)==0 display(['atom ' num2str(i)]);end
    
    whichindex =0;
    
    %in the case of the first voronoi type
    if i==1
        for j=1:numdifat
            vorhist(j).vorhist(1) = 0;
            vorhist(j).id(1).id = atom(i).vorind;
        end
    end
    
    %try to see if a pair matches existing vorind types
    for j=1:indexcount
        if atom(i).vorind == vorhist(1).id(j).id
            whichindex = j;
            break;
        end
    end
    
    %if it does not match, then make a new 
    if whichindex == 0
        indexcount = indexcount+1;
        for j=1:numdifat
            vorhist(j).vorhist(indexcount) = 0;
            vorhist(j).id(indexcount).id = atom(i).vorind;
        end
        whichindex = indexcount;
    end
    
    atomtype = atom(i).type;
    
    vorhist(atomtype).vorhist(whichindex) = vorhist(atomtype).vorhist(whichindex)+1;
    
end

for j=1:numdifat
    temphist = vorhist;
    maxhist = vorhist(j).vorhist;

    for i=1:length(maxhist)
        [val index]=max(maxhist);
        vorhist(j).vorhist(i) = temphist(j).vorhist(index);
        vorhist(j).id(i).id = temphist(j).id(index).id;
        maxhist(index) = -1;
    end
end
    
        