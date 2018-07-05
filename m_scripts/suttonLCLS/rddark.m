function [dark,images,imagnums]=rddark(scanno)

%%{
%read in an average dark
expno = 26410;
fina = runexpNO2fina(scanno,expno);
% scaninput.fina = fina;
% scaninput.DataSets = {'IPM2','IPM2TS'};
% scan = rdXPPdata(scaninput);
% exps=length(scan.scanvec);
dataset = '/Configure:0000/Run:0000/';
[scstepnames] = H5getobjnames(fina,dataset);
exps = length(scstepnames);

scanStepIndex = 1;
frameIndex = 1;

images=[];
imagnums=[];
for scanStepIndex = 1:exps
    try
        rawimage = double(rdPrincetonData(fina, scanStepIndex, frameIndex));
        images = cat(3, images, rawimage);
        imagnums = [imagnums scanStepIndex];
        fprintf( 'frame: %d\r',scanStepIndex);
    catch
       disp(['no frame: ' num2str(scanStepIndex)]);
    %fiend;
    end;
end; 
disp(['no of frames =' num2str(length(imagnums))]);
dark=mean(images,3);
end
%}