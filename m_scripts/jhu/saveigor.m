function saveigor(data, filename, header)
%SAVEIGOR.M
%
%Usage: saveigor(data, filename, {header});
%
%Saves a 1 or 2 dimentional column matrix 'data' to file 'filename'.
%The program generates an ASCII text file whose columns are exactly the
%same as those stored in the Matlab variable given to the function. If
%'header' is specified, its contents will be written to the first line.
%This header can be used to label each column with an Igor Pro wave name.
%
%Input: data - 1 or 2 dimentional matrix of values to be printed to the file
%       filename - string containing the name of the file to be created/
%               written.
%       header - {optional} single string containing space- or tab-separated 
%               column names that correspond to the columns in 'data'.
%               These will not necessarily line up above the column, but
%               Igor will still label the column waves appropriately.  
%
%Example: saveigor([time temp], 'temperature.txt', 'Time_axis Temperature')
%
%SH 11-14-07

if nargin<3 header =[]; end

fid = fopen(filename, 'w');

if ~isempty(header)
    fprintf(fid, [header '\n']);
end

for i=1:size(data,1)
    for j=1:size(data,2)
        fprintf(fid, '%20.10f  ',data(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);