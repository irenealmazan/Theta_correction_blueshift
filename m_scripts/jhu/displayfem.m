function displayfem(rootname, whichimages, cbar)

%SH 4/06 
%displayfem(rootname, whichimages, cbar)
%
%function for displaying up to nine fem images of the type femAA.im
%
%input
%  rootname - string which is the root of the fem files 
%  whichimages - array up to 9 numbers
%  cbar - 'y' or 'n' to include colorbar of each image
%output - up to 9 images displayed at once in a window

max = length(whichimages);
colormap(gray);
plotspot =1;
str2 = [rootname '.imax'];
imax = load(str2);
str3 = [rootname '.vk'];
vk = load(str3);

if length(whichimages)==1

    letter = whichimages-1;
    temp1= floor(letter/26);
    temp2=65+letter-temp1*26;

    str = [rootname 65+temp1 temp2 '.im']
    tempim = load(str);
    pcolor(imax,imax,tempim);
    
    %contourf(tempim);
    axis square;
    set(findobj('Type','surface'),'EdgeColor', 'none');
    set(findobj('Type','patch'),'EdgeColor', 'none');
    %title(['image 1 (' 65+temp1 temp2 ') ' num2str(vk(letter+1,1)) ' A^{-1}']);
  
    if(cbar=='y') colorbar;end

end

if length(whichimages)>1
    
    for i=1:max
        
        letter = whichimages(i)-1;
        temp1= floor(letter/26);
        temp2=65+letter-temp1*26;

        str = [rootname 65+temp1 temp2 '.im'];
        tempim = load(str);
        subplot(3,3,plotspot);
        pcolor(imax,imax,tempim);
        %contourf(tempim);
        axis square;
        set(findobj('Type','surface'),'EdgeColor', 'none');
        set(findobj('Type','patch'),'EdgeColor', 'none');
        title(['image ' num2str(i) ' (' 65+temp1 temp2 ') ' num2str(vk(letter+1,1)) ' A^{-1}']);
        plotspot = plotspot+1;
        if(cbar=='y') colorbar;end
 
    end

end






