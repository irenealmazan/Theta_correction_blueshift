function [image,images,imagnums]=rdimages(scanno,dark,darkname,plotflag)
% Read images
global fina;
expno = 26410;
fina = runexpNO2fina(scanno,expno);
dataset = '/Configure:0000/Run:0000/';
[scstepnames] = H5getobjnames(fina,dataset);
exps = length(scstepnames);
scanStepIndex = 1;
frameIndex = 1;
images=[];
imagnums=[];
for scanStepIndex = 1:exps
    try
        rawimage = double(rdPrincetonData(fina, scanStepIndex, frameIndex))-dark;
        images = cat(3, images, rawimage);
        imagnums = [imagnums scanStepIndex];
        fprintf('frame: %d\r',scanStepIndex');
        
    catch
        fprintf('no frame: %d\r',scanStepIndex');
    end;
end; 
disp(['no of frames =' num2str(length(imagnums))]);
image=mean(images,3);
if(plotflag)
    nfigure('Princeton');
    clf;imagesc(image,[0 100]);
    set(gca, 'ydir', 'normal');
    axis image;
    xlabel('tth (pixels)');
    ylabel('tth\_perp (pixels)');
    title(['Run ' num2str(scanno) ': using ' num2str(length(imagnums))  ' images, subtracted dark(' darkname ')']);
    drawnow;
    nfigure('averages');
    clf;subplot(2,1,1);
    plot(mean(image,2));
    xlabel('tth\_perp (pixels)');
    ylabel('I(adu)/frame');
    title(['Run ' num2str(scanno) ': using ' num2str(length(imagnums))  ' images, subtracted dark(run 31)']);

    subplot(2,1,2);
    plot(mean(image));
    xlabel('tth (pixels)');
    ylabel('I(adu)/frame');
    drawnow;
end


