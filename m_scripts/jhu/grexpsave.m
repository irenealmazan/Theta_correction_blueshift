
function grexpsave(grexp, filename)

fid = fopen(filename,'w');

fprintf(fid, '%4.0f\n', size(grexp,1));
fprintf(fid, '%10.8f\n', grexp(2,1));
fprintf(fid, '%12.8f\n', grexp(size(grexp,1),1));


for i=1:size(grexp,1)
    fprintf(fid,'%12.9f\t%12.9f\n', grexp(i,1), grexp(i,2));
end
st = fclose(fid);