function newarray = wrap1d(origarray, offset)

%SH 11-10-08
%shifts origarray by the value of offset imposing a 1d periodic boundary
%condition

origarray = origarray(:);
len = length(origarray);

if offset==0 newarray = origarray; end

if offset < 0 offset = len + offset; end

if offset > 0
    
    newarray(1:offset ,1) = origarray(len-offset+1 : len);
    newarray(offset+1:len, 1) = origarray (1:len-offset);
    
end
