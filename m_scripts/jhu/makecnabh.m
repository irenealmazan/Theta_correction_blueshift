function [cnabh cnabhall]= makecnabh(pairs, r, halfedge)

%SH commented 9-12-08
%
% This function compiles statistical distribution function information from
% a configuration based on the CNA indices of pairs of atoms stored in the
% 'pairs' structure (calculated using cna.m).  
%
% input:
%   pairs - structure calculated using cna.m
%   r - user defined r axis for the RDFs (in angstr)
%   halfedge - .5 * (natoms/atden)^(1/3).  half of the simulation box edge.
%
% output:
%   cnabh - structure containing the following fields
%       pairtype - cna pair type for which these statistics were tabulated
%       bh - total instances of this cna pair type as a function of
%       separation distance r.  
%       bh11 - instances of this cna pairing between two atoms of type 1
%       bh21 - instances of this cna pairing between atoms of type 1 and
%       atoms of type 2
%       bh22 - instances of this cna pairing between two atoms of type 2
%       bh... - continue for as many partial pairs as are present in the
%       simulation
%   cnabhall - matrix containing total CNA distribution colums sorted in
%   order of decreasing max value.  The same order is reflected in the
%   above 'cnabh' output structure.

pairidcnt=1;
rspacing = r(2)-r(1);

%for i=1:length(pairs)
for i=1:length(pairs)
    if(mod(i,4000)==0) display(['working on pair ' num2str(i) ' of ' num2str(length(pairs))]);
    end
    
    whichpairid =0;
    %display([num2str(pairs(i).cna) ' ' num2str(pairs(i).dist)]);
    
    %in the case of the first pair
    if i==1
        cnabh(pairidcnt).pairtype = pairs(i).cna;
        cnabh(pairidcnt).bh = zeros(length(r),1);
    
        %this part hard wired for pd ni p
        cnabh(pairidcnt).bh11 = zeros(length(r),1);
        cnabh(pairidcnt).bh21 = zeros(length(r),1);
        cnabh(pairidcnt).bh22 = zeros(length(r),1);
        cnabh(pairidcnt).bh31 = zeros(length(r),1);
        cnabh(pairidcnt).bh32 = zeros(length(r),1);
        cnabh(pairidcnt).bh33 = zeros(length(r),1);
    end
    
    %try to see if a pair matches existing cna types
    for j=1:pairidcnt
        %[pairs(i).cna cnabh(j).pairtype]
        if pairs(i).cna == cnabh(j).pairtype
            whichpairid = j;
        end
    end
    
    %if it does not match, then make a new pair id
    if whichpairid ==0
        pairidcnt=pairidcnt+1;
        cnabh(pairidcnt).pairtype = pairs(i).cna;
        cnabh(pairidcnt).bh = zeros(length(r),1);
        
        %this part hard wired for pd ni p
        cnabh(pairidcnt).bh11 = zeros(length(r),1);
        cnabh(pairidcnt).bh21 = zeros(length(r),1);
        cnabh(pairidcnt).bh22 = zeros(length(r),1);
        cnabh(pairidcnt).bh31 = zeros(length(r),1);
        cnabh(pairidcnt).bh32 = zeros(length(r),1);
        cnabh(pairidcnt).bh33 = zeros(length(r),1);
        
        whichpairid = pairidcnt;
    end
    
    %now put the pair distance in the right bin for this pair id
    %distribution
    dist = pairs(i).dist * halfedge;
    
    bin = round( dist / rspacing )+1;
    
    if bin<=length(r)
        cnabh(whichpairid).bh(bin) = cnabh(whichpairid).bh(bin)+1;
        
        %this is again hard wired for pd ni p
        %pairs(i).atom(1).type
        %pairs(i).atom(2).type
        
        if (pairs(i).atom(1).type == 1) & (pairs(i).atom(2).type == 1)
            
            cnabh(whichpairid).bh11(bin) = cnabh(whichpairid).bh11(bin)+1;
        end
        if pairs(i).atom(1).type == 2 & pairs(i).atom(2).type == 2
            cnabh(whichpairid).bh22(bin) = cnabh(whichpairid).bh22(bin)+1;
            %display('hit 2 2');
        end
        if pairs(i).atom(1).type == 3 & pairs(i).atom(2).type == 3
            cnabh(whichpairid).bh33(bin) = cnabh(whichpairid).bh33(bin)+1;
            %display('hit 3 3');
        end
        
        if (pairs(i).atom(1).type == 2 & pairs(i).atom(2).type == 1) | ... 
                (pairs(i).atom(1).type == 1 & pairs(i).atom(2).type == 2)
            cnabh(whichpairid).bh21(bin) = cnabh(whichpairid).bh21(bin)+1;
        end
        if (pairs(i).atom(1).type == 3 & pairs(i).atom(2).type == 1) | ... 
                (pairs(i).atom(1).type == 1 & pairs(i).atom(2).type == 3)
            cnabh(whichpairid).bh31(bin) = cnabh(whichpairid).bh31(bin)+1;
        end
        if (pairs(i).atom(1).type == 2 & pairs(i).atom(2).type == 3) | ... 
                (pairs(i).atom(1).type == 3 & pairs(i).atom(2).type == 2)
            cnabh(whichpairid).bh32(bin) = cnabh(whichpairid).bh32(bin)+1;
        end
        
    end
end

tempbh = cnabh;
for i=1:length(cnabh)
    maxbh(i)=max(cnabh(i).bh);
end

for i=1:length(cnabh)
    [val index]=max(maxbh);
    cnabh(i).bh = tempbh(index).bh;
    
    cnabh(i).bh11 = tempbh(index).bh11;
    cnabh(i).bh21 = tempbh(index).bh21;
    cnabh(i).bh22 = tempbh(index).bh22;
    cnabh(i).bh31 = tempbh(index).bh31;
    cnabh(i).bh32 = tempbh(index).bh32;
    cnabh(i).bh33 = tempbh(index).bh33;
    
    cnabh(i).pairtype = tempbh(index).pairtype;
    cnabhall(:,i)= tempbh(index).bh;
    maxbh(index)=-1;
end

