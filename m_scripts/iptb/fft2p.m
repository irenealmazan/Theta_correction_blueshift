function F=fft2p(array)

if(matlabpool('size')==0) matlabpool; end
numproc= matlabpool('size');

[numy numx]=size(array);

%do y dimension first
numperproc = floor(numy/numproc);
leftovers = mod(numy, numproc);

for i=1:numproc
    ymin= (i-1)*numperproc+1;
    ymax= i*numperproc;
    yrange(i,:)=[ymin ymax];
end

yrange

spmd
    Fsplit=fft(array(yrange(labindex,1):yrange(labindex,2),:),[],2);
end

F=[];
for i=1:numproc
    F=[F; Fsplit{i}];
end

spmd
    Fsplit=fft(F(:,yrange(labindex,1):yrange(labindex,2)),[],1);
end

F=[];
for i=1:numproc
    F=[F Fsplit{i}];
end
