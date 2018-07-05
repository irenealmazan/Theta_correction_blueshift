function [v,m]= checkvariance_novoyles(rootname, num, imortiff, whichrot)

if(imortiff(1)=='t')
    
    k=[.1:.05:1.3];
    for i=1:num
        file = [rootname num2str(k(i), '%3.2f') '.tif'];
        
        tempim = double(imread(file));
        pix = length(tempim);
   
        mean1 =0;
        isqterm = 0;
        var=0;

        for j=1:pix*pix
            mean1 = mean1 + tempim(j);
        end

        mean1 = mean1 / (pix*pix);
        %display(['mean is ' num2str(mean1) ' of image ' num2str(i)]);

        for j=1:pix*pix
            tempd = tempim(j);
            var = var+ (tempd-mean1) * (tempd-mean1);
        end

        var = var / (pix*pix -1);
        var = var /(mean1^2);
        
        v(i,1) = var;
        m(i,1)=mean1;

        immat(:,:,i,1) = tempim;
    
    end
    
end

if(imortiff(1)=='i')
    
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
        var=0;

        for j=1:pix*pix
            mean1 = mean1 + tempim(j);
        end

        mean1 = mean1 / (pix*pix);
        %display(['mean is ' num2str(mean1) ' of image ' num2str(i)]);

        for j=1:pix*pix
            tempd = tempim(j);
            var = var+ (tempd-mean1) * (tempd-mean1);
        end

        var = var / (pix*pix -1);

        v(i,1) = var/(mean1^2);
        m(i,1)=mean1;

        immat(:,:,i,1) = tempim;
    end
end