function datalayer = extractccdinfo(data, filepath, startimage, extract, ... 
    prenum, postnum)

if nargin<6
    prenum = 'image';
    postnum = '.hdf';
end

if(ndims(data) == 3)  %this is an area scan

    dims = size(data);
    numpts = dims(1) * dims(2);
    
    datalayer = zeros([dims(1) dims(2)]);
    
    counter =startimage;
    h=waitbarpos(0, 'loading ccd images');
    
    for i= dims(1):-1:1
        for j= 1: dims(2)
            
            ccd = hdfread([filepath 'image' num2str(counter,'%1.5d') '.hdf'], ...
                '/entry1/data/data', 'Index', {[1 1],[1 1], [1024 1024]});
            ccd = double(ccd(:,2:size(ccd,2)-1));
            
            if strcmp(extract.mode, 'maximum')
                datalayer(i,j) = max(max(ccd));
            elseif strcmp(extract.mode, 'minimum')
                datalayer(i,j) = min(min(ccd));
            elseif strcmp(extract.mode, 'mean')
                datalayer(i,j) = mean(mean(ccd));
            elseif strcmp(extract.mode, 'roi')
                xmin = extract.xrange(1);
                xmax = extract.xrange(2);
                ymin = extract.yrange(1);
                ymax = extract.yrange(2);
                
                datalayer(i,j) = mean(mean(ccd(ymin:ymax, xmin:xmax)));
            end
                
            
            
            counter = counter +1;            
             
            waitbar((counter-startimage)/ numpts);
            
        end
    end
    close(h);
end


if(ndims(data) ==2 ) %this is a line scan
    
    
end