function [threebody]= calcthreebody(pairs,atom, halfedge, r1,dr, raxis, thetaaxis)

r1 = r1/halfedge;
dr = dr/halfedge;
raxis = raxis/halfedge;

npairs = length(pairs)
ccount =1;

for i=1:npairs
   dist = pairs(i).dist;

   if (dist<r1+dr && dist>r1-dr)
       counterpairs(ccount) = pairs(i);
       ccount = ccount +1;
   end
    
end

ccount=ccount -1;

draxis=raxis(2)-raxis(1);
dtheta = thetaaxis(2) - thetaaxis(1);
numrbins = length(raxis);
numthetabins=length(thetaaxis);

threebody(numrbins, numthetabins) = 0;

for i=1:length(atom)
    origincentr = atom(i).pos;
    if(mod(i,20)==0) display(['doing origin ' num2str(i) ' of ' num2str(length(atom))]);end
    
    for j=1:ccount
        countercentr = counterpairs(j).centroid;
        pairdist = sqrt((origincentr(1)-countercentr(1))^2 + ...
            (origincentr(2)-countercentr(2))^2 + ...
            (origincentr(3)-countercentr(3))^2);
        rplace = 1+round(pairdist/draxis);
    
        v1 = origincentr - countercentr;
        
        v2 = counterpairs(j).atom(1).pos - counterpairs(j).atom(2).pos;
        pairtheta = abs(dot(v1,v2)/(norm(v1)*norm(v2)));
        pairtheta = acosd(pairtheta);
        thetaplace = 1+round(pairtheta/dtheta);
        
        if(thetaplace<=numthetabins && rplace<=numrbins)
            threebody(rplace,thetaplace) = threebody(rplace, thetaplace)+1;
        end
        %pause;    
    end
end

threebody(1,1)=0;