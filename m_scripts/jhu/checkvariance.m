 function [v,m]=checkvariance(rootname, num, orient, imortiff, whichrot)


if( ~isa(rootname,'char'))
    for i=1:num

        tempim = rootname(:,:,i);

        pix = length(tempim);

        mean1 =0;
        isqterm = 0;

        for j=1:pix*pix
            tempim(j);
            mean1 = mean1 + tempim(j);
            isqterm = isqterm + tempim(j)*tempim(j);
        end

        mean1 = mean1 / (pix*pix);
        
        v(i,1) = (isqterm/(pix*pix)) / (mean1*mean1);
        v(i,1) = v(i,1)-1;
        m(i,1)=mean1;

        immat(:,:,i,1) = tempim;
    end
end 
 
if(imortiff(1) =='t' && isa(rootname,'char'))
    
    k=[.1:.05:1.3];
    for i=1:num
        file = [rootname num2str(k(i), '%3.2f') '.tif'];
        
        tempim = double(imread(file));
        pix = length(tempim);
   
        mean1 =0;
        isqterm = 0;

        for j=1:pix*pix
            tempim(j);
            mean1 = mean1 + tempim(j);
            isqterm = isqterm + tempim(j)*tempim(j);
        end

        mean1 = mean1 / (pix*pix);

        v(i,1) = (isqterm/(pix*pix)) / (mean1*mean1);
        v(i,1) = v(i,1)-1;
        m(i,1)=mean1;

        immat(:,:,i,1) = tempim;        
    
    end
    
end



if(imortiff(1)== 'i'&& isa(rootname,'char'))
    for i=1:num

        letter = i-1;
        temp1= floor(letter/26);
        temp2=65+letter-temp1*26;

        str = [rootname 65+temp1 temp2 '.im'];
        tempim = load(str);
    
        pix = size(tempim,2);
        
        tempim = tempim(1+pix*whichrot:whichrot*pix+pix, :);
        
        mean1 =0;
        isqterm = 0;

        for j=1:pix*pix
            tempim(j);
            mean1 = mean1 + tempim(j);
            isqterm = isqterm + tempim(j)*tempim(j);
        end

        mean1 = mean1 / (pix*pix);
        
        v(i,1) = (isqterm/(pix*pix)) / (mean1*mean1);
        v(i,1) = v(i,1)-1;
        m(i,1)=mean1;

        immat(:,:,i,1) = tempim;
    end

    if(orient>1)
        for i=1:num

            letter = i-1;
            temp1= floor(letter/26);
            temp2=65+letter-temp1*26;

            str = [rootname 65+temp1 temp2 '.im2'];
            tempim = load(str);

            pix = length(tempim);

            mean1 =0;
            isqterm = 0;

            for j=1:pix*pix
                tempim(j);
                mean1 = mean1 + tempim(j);
                isqterm = isqterm + tempim(j)*tempim(j);
            end

            mean1 = mean1 / (pix*pix);

            v(i,2) = (isqterm/(pix*pix)) / (mean1*mean1);
            v(i,2) = v(i,2)-1;
            m(i,2)=mean1;

            immat(:,:,i,2) = tempim;
        end
    end

    if(orient>2)    
        for i=1:num

            letter = i-1;
            temp1= floor(letter/26);
            temp2=65+letter-temp1*26;

            str = [rootname 65+temp1 temp2 '.im3'];
            tempim = load(str);

            pix = length(tempim);

            mean1 =0;
            isqterm = 0;

            for j=1:pix*pix
                tempim(j);
                mean1 = mean1 + tempim(j);
                isqterm = isqterm + tempim(j)*tempim(j);
            end

            mean1 = mean1 / (pix*pix);

            v(i,3) = (isqterm/(pix*pix)) / (mean1*mean1);
            v(i,3) = v(i,3)-1;
            m(i,3)=mean1;

            immat(:,:,i,3) = tempim;
        end
    end

    if(orient>3)
        for i=1:num

            letter = i-1;
            temp1= floor(letter/26);
            temp2=65+letter-temp1*26;

            str = [rootname 65+temp1 temp2 '.im4'];
            tempim = load(str);

            pix = length(tempim);

            mean1 =0;
            isqterm = 0;

            for j=1:pix*pix
                tempim(j);
                mean1 = mean1 + tempim(j);
                isqterm = isqterm + tempim(j)*tempim(j);
            end

            mean1 = mean1 / (pix*pix);

            v(i,4) = (isqterm/(pix*pix)) / (mean1*mean1);
            v(i,4) = v(i,4)-1;
            m(i,4)=mean1;

            immat(:,:,i,4) = tempim;
        end
    end

    v(:,5)= sum(v(:,1:orient),2)/orient;
    m(:,5)= sum(m(:,1:orient),2)/orient;
end