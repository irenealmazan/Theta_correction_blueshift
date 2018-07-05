function out = mirrorobject(in, dims);

if nargin<2 dims=[1 2 3]; end
    
tempmat = in;

if(ismember(1,dims))
    
    s=size(in,1);
    for i=1:s
        temp = tempmat(i,:,:);
        out((s+1)-i,:,:)=temp;
    end
    tempmat = out;
 
end

if(ismember(2,dims))

    s=size(in,2);
    for i=1:s
        temp = tempmat(:,i,:);
        out(:,(s+1)-i,:)=temp;
    end
    tempmat = out;
end

if(ismember(3,dims))
    
    s=size(in,3);
    for i=1:s
        temp = tempmat(:,:,i);
        out(:,:,(s+1)-i)=temp;
    end
end
