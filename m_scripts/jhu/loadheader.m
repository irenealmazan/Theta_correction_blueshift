function array = loadheader(filename, skiplines)

fid1=fopen(filename, 'r');

%skip the header lines
for i=1:skiplines
    temp = fgetl(fid1);
end

%find out how many columns there are
temp = str2num(fgetl(fid1));
columns = length(temp);

fs=[];
for i=1:columns
    fs=[fs '%f '];
end

array = fscanf(fid1, fs, [columns inf]);
array = array';

fclose(fid1);

array( 2:size(array,1)+1 ,:) = array;
array( 1,:) = temp;