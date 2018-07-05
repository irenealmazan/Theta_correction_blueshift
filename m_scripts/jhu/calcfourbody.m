function [fourbody]= calcfourbody(pairs, halfedge, r1,r3,dr, raxis, thetaaxis)

r1 = r1/halfedge;
r3 = r3/halfedge;
dr = dr/halfedge;
raxis = raxis/halfedge;

npairs = length(pairs)
ocount =1;
ccount =1;

for i=1:npairs
   dist = pairs(i).dist;
   if (dist<r1+dr & dist>r1-dr) 
       originpairs(ocount) = pairs(i);
       ocount = ocount +1;
   end

   if (dist<r3+dr & dist>r3-dr)
       counterpairs(ccount) = pairs(i);
       ccount = ccount +1;
   end
    
end

ocount=ocount -1;
ccount=ccount -1;

draxis=raxis(2)-raxis(1);
dtheta = thetaaxis(2) - thetaaxis(1);
numrbins = length(raxis);
numthetabins=length(thetaaxis);

fourbody(numrbins, numthetabins) = 0;
h = waitbar(0,'calculating four body correlation');

for i=1:ocount
    origincentr = originpairs(i).centroid;
    v1 = originpairs(i).atom(1).pos(1:3) - originpairs(i).atom(2).pos(1:3);
    normv1 = norm(v1);
    %if(mod(i,20)==0) display(['doing origin ' num2str(i) ' of ' num2str(ocount)]);end
    
    for j=1:ccount
        countercentr = counterpairs(j).centroid;
        pairdist = sqrt((origincentr(1)-countercentr(1))^2 + ...
            (origincentr(2)-countercentr(2))^2 + ...
            (origincentr(3)-countercentr(3))^2);
        rplace = 1+round(pairdist/draxis);
       
        v2 = counterpairs(j).atom(1).pos(1:3) - counterpairs(j).atom(2).pos(1:3);
        pairtheta = dot(v1,v2)/(normv1*norm(v2));
        pairtheta = acosd(pairtheta);
        thetaplace = 1+round(pairtheta/dtheta);
        
        if(thetaplace<=numthetabins && rplace<=numrbins)
            fourbody(rplace,thetaplace) = fourbody(rplace, thetaplace)+1;
        end
        %pause;    
    end
    waitbar(i/ocount)
end

close(h);

fourbody(1,1)=0;
imagesc(fourbody);axis equal;