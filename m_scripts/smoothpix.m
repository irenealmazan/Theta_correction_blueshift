function smoothccd=smoothpix(ccd, hotpixels)

%SH 4-23-09
%
%This function takes a pixel index designated by the hotpixels matrix and
%re-assigns its value as a mean of its 8 neighboring pixels

smoothccd = ccd;
col = [1:size(ccd,1)];
row = [1:size(ccd,2)];

for i = 1:length(hotpixels)
    
    if(numel(hotpixels(i,:))==1)
        [pixi pixj] = ind2sub(size(ccd), hotpixels(i));
        pix=[pixi pixj];
    else
        pix = hotpixels(i,:);
    end
    
    %[i hotpixels(i) pix]
    %imagesc(ccd(pix(1)-2:pix(1)+2, pix(2)-2:pix(2)+2));colorbar;pause
    
    
    ave=[];
    %left of center column
    if (ismember(pix(1)-1,col) && ismember(pix(2)-1, row)) 
        ave = [ave ccd(pix(1)-1, pix(2)-1)]; end
    if (ismember(pix(1),col) && ismember(pix(2)-1, row)) 
        ave = [ave ccd(pix(1), pix(2)-1)]; end
    if (ismember(pix(1)+1,col) && ismember(pix(2)-1, row)) 
        ave = [ave ccd(pix(1)+1, pix(2)-1)]; end
    %right of center column
    if (ismember(pix(1)-1,col) && ismember(pix(2)+1, row)) 
        ave = [ave ccd(pix(1)-1, pix(2)+1)]; end
    if (ismember(pix(1),col) && ismember(pix(2)+1, row)) 
        ave = [ave ccd(pix(1), pix(2)+1)]; end
    if (ismember(pix(1)+1,col) && ismember(pix(2)+1, row)) 
        ave = [ave ccd(pix(1)+1, pix(2)+1)]; end
    %above and below the pixel
    if (ismember(pix(1)-1,col) && ismember(pix(2)+1, row)) 
        ave = [ave ccd(pix(1)-1, pix(2))]; end
    if (ismember(pix(1)+1,col) && ismember(pix(2)+1, row)) 
        ave = [ave ccd(pix(1)+1, pix(2))]; end
    
    %[ccd(pix(1),pix(2)) ccd(hotpixels(i))]
    %ave
    
    %if(length(ave)==8)
    %    imagesc(ccd(pix(1)-2:pix(1)+2, pix(2)-2:pix(2)+2));colorbar;pause
    %end
    
    ave=double(ave);
    sdev = std(ave);
    ave(find(ave>mean(ave)+sdev)) = [];
    ave(find(ave<mean(ave)-sdev)) = [];
    
    ave = mean(ave);
    smoothccd(pix(1), pix(2))=ave;
        
%     if(pix(1)==203 && (pix(2)==369)) 
%         imagesc(ccd(pix(1)-2:pix(1)+2, pix(2)-2:pix(2)+2));colorbar;
%         ave=ccd(pix(1)-1:pix(1)+1, pix(2)-1)';
%         ave = [ave ccd(pix(1)-1:pix(1)+1, pix(2)+1)'];
%         ave = [ave ccd(pix(1)-1, pix(2)) ccd(pix(1)+1, pix(2))];
%         ave
%         pause; 
%     end
end

