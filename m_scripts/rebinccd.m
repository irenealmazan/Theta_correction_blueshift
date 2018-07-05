function binnedimage = rebinccd(oldimage, binrate, flag)

%SH 11-7-08
%
%can handle 2 or 3 dimensional arrays.  if 3d, then the 3rd axis is
%unbinned and the binning is just in the x - y plane.  in the case where
%the x or y pixel axes are indivisble by binrate, then the remainder will
%be evened out by eliminating the first one (or more) columns or rows as
%appropriate. to further reduce size, the output array are single precision
%floats rather than double
%
%input oldimage - 2 or 3 dimensional array to be reduced
%       binrate - the size of the 2d square of pixels that will be comined
%       into one.  e.g. if binrate = 2 then a 2x2 area of the old image
%       will be averaged into one pixel in the new image.  
%
%output binnedimage - the result of the program 

%if not indicated, use a sum of the pixel vals rather than mean
if(nargin<3) flag='s'; end

[s1 s2 s3] = size(oldimage);

extrarows = mod(s1,binrate);
extracols = mod(s2, binrate);

if extrarows~=0
    oldimage(1:extrarows,:,:)= [];
    display(['removed ' num2str(extrarows) ' rows from array']);
end

if extracols ~=0
    oldimage(:, 1:extracols,:) = [];
    display(['removed ' num2str(extracols) ' columns from array']);
end

[s1 s2 s3] = size(oldimage);
col = 0;
row =0;
element = zeros(1, binrate*binrate);
binnedimage = zeros(s1/binrate, s2/binrate, s3);

h=waitbarpos(0, 'binning individual CCD images');

for s3dir = 1:s3
    
    col=0;
    for i=1:binrate:s2 %go down the columns

        col = col+1;
        %display(['column ' num2str(col)]);
        row =0;
        for j=1:binrate:s1 % go down the rows

            row = row+1;
            %display(['row ' num2str(row)]);

            counter = 1;
            for k=0:binrate-1 %defines the colums to be binned
                for l=0:binrate-1 %defines the rows to be binned

                    element(counter) = oldimage(j+l ,i+k ,s3dir);
                    if element(counter) == nan element(counter)=0; end
                    if element(counter) == Inf element(counter)=0; end
                    if element(counter) == -Inf element(counter)=0; end

                    %display(['adding to element ' num2str(j+l) '  ' num2str(i+k)]);

                    counter = counter +1;
                end
            end

            %element
            if(strcmp('m',flag(1)))
                binel = mean(element);
            else
                binel = sum(element);
            end
            
            %display (['add to row ' num2str(row) '  and column ' num2str(col)]);
            binnedimage(row,col,s3dir) = single(binel);
           %pause 
        end
    end
    waitbar(s3dir/s3);
end

close(h);