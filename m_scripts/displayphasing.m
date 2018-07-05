function displayphasing(A,B,C,D, support, dim, hfig, hchi)

figure(hfig);
clf 


if dim==3
    
    x=[1:size(A,2)];
    y=[1:size(A,1)];
    z=[1:size(A,3)];

    [X,Y,Z]= meshgrid(x,y,z);

    subplot(4,2,1);
    displayisosurf( fftshift((abs(A))), -.02, 'g');
    subplot(4,2,2);
    imagesc(fftshift(angle(A(:,:,round(size(A,3)/2)))));axis equal;

    subplot(4,2,3);
    displayisosurf( fftshift((abs(B))), -.02, 'g');
    %displayisosurf( fftshift(log(abs(B))), -.95, 'g');
    subplot(4,2,4);
    imagesc(fftshift(angle(B(:,:,round(size(A,3)/2)))));axis equal;

    subplot(4,2,5);
    displayisosurf( support, .999, 'y');
    alpha(.3);
    displayisosurf( abs(C), -.1, 'g');
    subplot(4,2,6);
    imagesc(angle(C(:,:,round(size(A,3)/2))));axis equal;

    subplot(4,2,7);
    displayisosurf( abs(D), -.5, 'g');
    subplot(4,2,8);
    imagesc(angle(D(:,:,round(size(A,3)/2))));axis equal;

    pause(.01)
       
end

if dim ==2
    for i=1:2:size(A,3)

        subplot(4,2,1);
        imagesc(fftshift(log(abs(A(:,:,i)))));
        subplot(4,2,2);
        imagesc(fftshift(angle(A(:,:,i))));

        subplot(4,2,3);
        imagesc(fftshift(log(abs(B(:,:,i)))));
        subplot(4,2,4);
        imagesc(fftshift(angle(B(:,:,i))));

        subplot(4,2,5);
        imagesc(abs(C(:,:,i)));
        subplot(4,2,6);
        imagesc(angle(C(:,:,i)));

        subplot(4,2,7);
        imagesc(abs(D(:,:,i)));
        subplot(4,2,8);
        imagesc(angle(D(:,:,i)));

        pause(.01)
    end
end

figure(hchi);