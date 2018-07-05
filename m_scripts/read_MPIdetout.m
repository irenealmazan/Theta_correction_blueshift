function  [MPI_I MPIHmat MPIKmat MPILmat] = read_MPIdetout(filename)

fid = fopen(filename,'r');

line = fgetl(fid);% blank
line = fgetl(fid); %save date and time
line = fgetl(fid); %contains the SPEC command
remain = line;
nscansteps='junk';
while remain
    [nscansteps remain] = strtok(remain);
    %display(nscansteps);
    if strcmp(nscansteps,'of') break; end
end
[nscansteps remain] = strtok(remain);
nscansteps = str2num(nscansteps);

line = fgetl(fid); %motor positions
line = fgetl(fid); %evaluating sample from file
remain=line;
for ii=1:5
    [sampfile remain] = strtok(remain);
end
line = fgetl(fid); %contains detector info
remain=line;
for ii=1:26
    [token remain] = strtok(remain);
    %display([num2str(ii) '  ' token]);
    if (ii==6) Hcen=(token); end
    if (ii==7) Kcen=(token); end
    if (ii==8) Lcen=(token); end
    if (ii==11) xpix=str2num(token); end
    if (ii==15) ypix=str2num(token); end
    if (ii==19) pxwidth=str2num(token); end
    if (ii==26) detdist=str2num(token); end
end
fclose(fid);

MPI_I=zeros(ypix, xpix, nscansteps);
MPIHmat=zeros(ypix,xpix, nscansteps);
MPIKmat=zeros(ypix,xpix, nscansteps);
MPILmat=zeros(ypix,xpix, nscansteps);
MPIqxmat=zeros(ypix,xpix, nscansteps);
MPIqymat=zeros(ypix,xpix, nscansteps);
MPIqzmat=zeros(ypix,xpix, nscansteps);
counter=1;
i=sqrt(-1);

fid = fopen(filename,'r');

for scan= 1:nscansteps
    
    for ii=1:6 line=fgetl(fid); end %skip the first 6 lines of info
    
    %now get the detector info
    clear mat
    mat = single(fscanf(fid, '%g %g %g %g %g %g %g %g', [8 xpix*ypix]));
    
    abstemp = mat(4,:); abstemp= reshape(abstemp, xpix, ypix); abstemp = abstemp';
    angtemp = mat(5,:); angtemp= reshape(angtemp, xpix, ypix); angtemp = angtemp';
    MPI_I(:,:,scan) = abstemp .* exp(i*angtemp);
    
    temp = mat(6,:); temp= reshape(temp, xpix, ypix); 
    MPIHmat(:,:,scan) = temp';
    temp = mat(7,:); temp= reshape(temp, xpix, ypix); 
    MPIKmat(:,:,scan) = temp';
    temp = mat(8,:); temp= reshape(temp, xpix, ypix); 
    MPILmat(:,:,scan) = temp';
    
    temp = mat(1,:); temp= reshape(temp, xpix, ypix); 
    MPIqxmat(:,:,scan) = temp';
    temp = mat(2,:); temp= reshape(temp, xpix, ypix); 
    MPIqymat(:,:,scan) = temp';
    temp = mat(3,:); temp= reshape(temp, xpix, ypix); 
    MPIqzmat(:,:,scan) = temp';
    line=fgetl(fid);
    
    
end

mid1 = ceil((size(MPI_I,1)+1)/2);
mid2 = ceil((size(MPI_I,2)+1)/2);
mid3 = 1;  if(size(MPI_I,3)>1) mid3=ceil((size(MPI_I,3)+1)/2); end

MPIH=MPIHmat(:,mid2,mid3);
MPIL=MPILmat(mid1,:,mid3);
clear mid1 mid2 mid3

MPIHmat=squeeze(MPIHmat);
MPIKmat=squeeze(MPIKmat);
MPILmat=squeeze(MPILmat);
MPI_I=squeeze(MPI_I);
