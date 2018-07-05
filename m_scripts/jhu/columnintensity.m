function columnintensity

inten =0;
diam=9.38;
Q=.61/(diam/2);
k=[.2:.02:1.2]';
r = [0:.1:diam+1]';
bh(length(r),1)=0;
bhwt(length(r),1)=0;  

rprimespacing =.01;
rprime=[0:rprimespacing:10]'; 
gr = zeros(length(rprime),1);
grnss = zeros(length(rprime),1);

pos=rand(150,2)*diam-diam/2;
posmag=sqrt(pos(:,1).^2 + pos(:,2).^2);

pos=pos(find(posmag<=diam/2),:);
posmag = posmag(find(posmag<=diam/2));

plot(pos(:,1),pos(:,2),'o');axis equal;

Igr = zeros(length(k),1);
Igrnss = zeros(length(k),1);
Ihcnff = zeros(length(k),1);
Ihcnffnss = zeros(length(k),1);
Ihcgr = zeros(length(k),1);
Idf =zeros(length(k),1);
Ixray = zeros(length(k),1);
Ixraynss = zeros(length(k),1);

for i=1:length(pos) %use each origin
    sigi = posmag(i);
    appfunci = besselj(1,2*pi*Q*sigi)/(2*pi*Q*sigi);
    
    for j=1:length(pos)
        sigj = posmag(j);
        appfuncj = besselj(1,2*pi*Q*sigj)/(2*pi*Q*sigj); 
        rij = sqrt((pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2);
        
        Ihcnff = Ihcnff + appfunci*appfuncj * besselj(0,2*pi*k*rij);
        
        binadd = round(rij/rprimespacing)+1;
        gr(binadd) = gr(binadd)+ appfunci*appfuncj;
        
        if(rij>0) 
            Ihcnffnss = Ihcnffnss + appfunci*appfuncj * besselj(0,2*pi*k*rij); 
            Ixraynss = Ixraynss + sin(k*rij*2*pi)./(k*2*pi*rij);
            Ixray = Ixray+1;
            grnss(binadd) = grnss(binadd)+ appfunci*appfuncj;
        end

    end
end

for i=1:length(k)
    Igr(i) = sum(gr .* besselj(0,2*pi*k(i)*rprime));
    Igrnss(i) = sum(grnss .* besselj(0,2*pi*k(i)*rprime));
end

plot(k,Ihcnff, k,Ihcnffnss,k,Igr,'--', k,Igrnss,'--');
%figure
%plot(k,Ixraynss,'--');
%figure
%plot(rprime,gr, rprime, grnss);
%figure 
%plot(k,Igr, k,Igrnss);

% for i=1:1
%     
%     sigi = sqrt((pos(i,1)*halfedge-pos(i,2)*halfedge)^2 + (imx-imy)^2 ); 
%     
%     appfunci = (Q/sigi) * besselj(1, 2*pi*sigi*Q);
%         
%     for j=1:length(pos)
%         rij = sqrt((pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2 + ...
%             (pos(i,3)-pos(j,3))^2);
%         rij = rij*halfedge;
%                 
%         sigij = sqrt((pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2 );
%         sigij = sigij*halfedge;
%         
%         sigj = sqrt((pos(j,1)*halfedge-pos(j,2)*halfedge)^2 + (imx-imy)^2 ); 
%     
%         appfuncj = (Q/sigj) * besselj(1, 2*pi*sigj*Q);
% 
%         bh(round(rij/.1)+1) = bh(round(rij/.1)+1) + 1;
%         bhwt(round(rij/.1)+1) = bhwt(round(rij/.1)+1) + appfuncj * appfunci;        
%         
%         inten = inten+ appfunci * appfuncj * besselj(0, 2*pi*k*sigij);
%         
%     end
% end
% 
% plot(r,bh, r, bhwt);