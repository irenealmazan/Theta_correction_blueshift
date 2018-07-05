function out = zeropad3d(in, outcol, outrow, outheight, coloffset, rowoffset, heightoffset)

%SH 8-3-09
%
%this function pads a 3-d array with zeros, while the original array
%remains centered.  this is not how fft2 handles zero padding, which works
%by simply adding zeros in +x and +y only.

if nargin<5
    coloffset=0;
    rowoffset=0;
    heightoffset=0;
end

if (size(in, 2) > outcol)
    flag=1;
    while flag
        in(:,1,:)=[];
        if size(in,2) <= outcol break; end
        in(:,end,:) = [];
        if size(in,2) <= outcol break; end
%         if size(in,2) <= outcol
%             flag=0; 
%             display('shrink number of columns');
%         end
    end
end

if (size(in, 1) > outrow)
    flag=1;
    while flag
        in(1,:,:)=[];
        if size(in,1) <=outrow break; end
        in(end,:,:) = [];
        if size(in,1) <=outrow break; end
%         if size(in,1) <= outrow
%             flag=0; 
%             display('shrink number of rows');
%         end
    end
end

if (size(in, 3) > outheight)
    flag=1;
    while flag
        in(:,:,1)=[];
        if size(in,3) <=outheight break; end
        in(:,:,end) = [];
        if size(in,3) <=outheight break; end
        if size(in,3) <= outheight
            flag=0; 
            display('shrink height');
        end
    end
end

[inrow incol inheight] = size(in);

out(outrow, outcol, outheight)=0;

lowcol = outcol/2 - incol/2 +1;
highcol = outcol/2 + incol/2;

lowrow = outrow/2 - inrow/2 +1;
highrow = outrow/2 + inrow/2;

lowheight = outheight/2 - inheight/2+1;
highheight = outheight/2 + inheight/2;

%display(num2str(size(in)));
%display(num2str([lowrow highrow lowcol highcol lowheight highheight]));
out([floor(lowrow):floor(highrow)] + rowoffset, [floor(lowcol):floor(highcol)] + coloffset, [floor(lowheight):floor(highheight)] + heightoffset) = in;
