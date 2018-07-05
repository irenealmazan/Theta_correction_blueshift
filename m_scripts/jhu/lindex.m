function [L, tree] = lindex(atom, cnaid)


for l=1:length(cnaid) %start a tree from each 

    %make the head node of this tree----------------------------------
    
    clear node;  %get rid of any old tree
    node(1).id = cnaid(l);  %atomic id
    node(1).generation = 0;  %indicates that this is the head of the tree
    node(1).parent = []; %again, its a head
    node(1).upbonds = []; %same

    %now fill in the fields that require calculatoin
    nn = atom(node(1).id).nn;  %take the nn from main structure
    
    clear nnlisttemp
    
    nnlisttemp = atom(node(1).id).neighborlist;
    
    node(1).neighbors = intersect(nnlisttemp, cnaid); %limits it to neighbors in the cna analysis
    node(1).nn = length(node(1).neighbors);
    
    if ~isempty(node(1).neighbors)
        for m=1:node(1).nn %this is all possible bonds to relevant nn
            node(1).nnbonds(m) = sqrt(node(1).id) * sqrt(node(1).neighbors(m));
        end
    else
        node(1).nnbonds = [];
    end
    
    %in the case of the head, the allowed offspring are all
    %the possible bonds.  no restrictions on backtracking
    %yet, allowing us to do:
    node(1).downbonds = node(1).nnbonds;
    node(1).offspring = node(1).neighbors;

    
    
    
    %make the next gerations as needed  ---------------------------------
    
    growing =1;
    generation =1; %the generation we are on
    nodecnt =2;
    
    while growing
        
        %find the parents that can potentially spawn into this generation
        clear parents
        cnt=1;
        for m=1:nodecnt-1  %scan all existing nodes
            gen = node(m).generation;
            if gen == generation -1
                parents(cnt) = m;
                cnt=cnt+1;
            end
        end
        
        thisgen=[];
        thisgencnt=0;
        
        for m=1:length(parents) %create this parent's branch
            
            kids = node(parents(m)).offspring;
            
            if ~isempty( node(parents(m)).offspring )
            
                for n=1:length(kids)
                    node(nodecnt).id = kids(n);
                    node(nodecnt).generation = generation;
                    node(nodecnt).parent = parents(m);

                    %figure out the upbonds, aka bonds following up to the node
                    thisupbond = sqrt(node(nodecnt).id) * sqrt(node(parents(m)).id);
                    node(nodecnt).upbonds = cat(2,node(parents(m)).upbonds, thisupbond);


                    %first try to see if the nn work has been done
                    match=[];
                    for o=1:length(nodecnt-1) %check all existing nodes
                        %if ~isempty(find(node(o).id == node(nodecnt).id))
                        if node(o).id == node(nodecnt).id
                            match = o;
                            break
                        end
                    end

                    if ~isempty(match)
                        node(nodecnt).nn = node(match).nn;
                        node(nodecnt).neighbors = node(match).neighbors;
                        node(nodecnt).nnbonds = node(match).nnbonds;
                    end


                    %figure out the neighbors the long way if we must

                    if isempty(match) %no suitable matches were found

                        clear nnlisttemp
                        nn = atom(node(nodecnt).id).nn;  %take the nn from main structure
                        
                        nnlisttemp = atom(node(nodecnt).id).neighborlist;

                        node(nodecnt).neighbors = intersect(nnlisttemp, cnaid); %limits it to neighbors in the cna analysis
                        node(nodecnt).nn = length(node(nodecnt).neighbors);

                        for o=1:node(nodecnt).nn %this is all possible bonds to relevant nn
                            node(nodecnt).nnbonds(o) = sqrt(node(nodecnt).id) * ... 
                                sqrt(node(nodecnt).neighbors(o));
                        end
                    end                

                    
                    %kids gotta have rules. we will pick out what offspring are
                    %allowed according to two criteria:
                    %1. you cant produce an offspring that is the same identity
                    %as your parent
                    %2. once all your potential downbonds have been used in
                    %your upbonds history, you are castrated and cant have more
                    %offspring.  aka this branch stops growing

                    %rule number 1
                    parid = node(parents(m)).id;
                    offspring = node(nodecnt).neighbors(find(node(nodecnt).neighbors~=parid));
                    
                    %offspring
                    
                    %rule number 2
                    downbonds = sqrt(offspring) .* sqrt(node(nodecnt).id);

                    nono = intersect(downbonds, node(nodecnt).upbonds);
                    
                    %downbonds
                    %nono
                    %node(nodecnt)
                    
                    if ~isempty(nono)
                        for o=1:length(nono)
                            offspring = offspring(find(downbonds ~= nono(o)));
                            downbonds = downbonds(find(downbonds ~= nono(o)));
                        end
                    end

                    node(nodecnt).offspring = offspring;
                    node(nodecnt).downbonds = downbonds;

                    thisgencnt = thisgencnt +1;
                    thisgen(thisgencnt)=nodecnt;
                    nodecnt = nodecnt+1;
                end
            
            end %end if ~isempty(offspring)
            
        end
    
        temp=[];
        for m=1:length(thisgen)
             temp= cat(2,temp,node(thisgen(m)).offspring);
        end
        if isempty(temp)
            growing=0;
        end
        
%         thisgen
%         if generation ==5
%             growing =0;
%         end
        generation = generation+1;
    end %growing
    L(l) = generation -1;
    tree(l).node = node;
end
