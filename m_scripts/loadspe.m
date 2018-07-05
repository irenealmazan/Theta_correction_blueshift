function [flag,data,xval,yval] = loadspe(filepathname)
fid=fopen(filepathname,'r');
if fid == 0
    errordlg('Can''t open this SPE file.','Error','modal');
    flag = -1;
    data = [];
    return;
else
    try
        fseek(fid,42,'bof');
        xdim = fread(fid,1,'uint16');
        fseek(fid,656,'bof');
        ydim = fread(fid,1,'uint16');
        fseek(fid,1446,'bof');
        nframes = fread(fid,1,'int32');
        fseek(fid,108,'bof');
        datatype = fread(fid,1,'int16');
        fseek(fid,362,'bof');
        xval = fread(fid,8,'char')';
        xval = char(xval);
        xval = str2num(xval);
        fseek(fid,373,'bof');
        yval = fread(fid,9,'char')';
        yval = char(yval);
        yval = str2num(yval);
        
%         header=fread(fid,2050,'uint16');% 4100 bytes/2
        fseek(fid,4100,'bof');
        switch datatype
            case 0
                ImMat=fread(fid,'*float32');
            case 1
                ImMat=fread(fid,'*int32');
            case 2
                ImMat=fread(fid,'*int16');
            case 3
                ImMat=fread(fid,'*uint16');
        end
        data=reshape(ImMat,xdim,ydim,nframes);
        fclose(fid);
        flag = 1;
    catch
        errordlg('Can''t open this SPE file.','Error','modal');
        flag = -1;
        data = [];
        return;
    end
end
