function [energy,spectra]=rdhdf(name,dis)
% [energy,spectra]=rdhdf(name,dis)
% name is filename; can be super style (fred.103) or
% hdf file (fred103.hdf).  dispara is optional, if missing, global
% parameters are displayed
% read in a super .hdf file and output energy and mca trace
% 17mar06SMB

dispara=0;
if nargin<2,dispara=1;end

if isempty(strfind(name,'.hdf')),
    idx=findstr(name,'.')
    hdfname=[name(1:idx-1) name(idx+1:idx+3) '.hdf'];
else
    hdfname=name;
end
if(dispara),disp(hdfname);end

sd_id=hdfsd('start',hdfname,'read');
[nsets,nglob]=hdfsd('fileinfo',sd_id);

startval=0;endval=0;eslope=0;eint=0;

for jj=1:nglob
    paraname = hdfsd('attrinfo',sd_id,jj-1);
    para= hdfsd('readattr',sd_id,jj-1);
    if(dispara),disp([paraname ' = ' num2str(para)]);end
    if strcmp(paraname,'mca_start'),startval=double(para);end
    if strcmp(paraname,'mca_end'),endval=double(para);end
    if strcmp(paraname,'energy_xscale'),eslope=double(para);end
    if strcmp(paraname,'energy_yint'),eint=double(para);end
end

range=endval-startval+1;

spectra=zeros(range,nsets);
energy=eslope*[startval:endval]+eint;

if(dispara),disp(nsets);end
start=[0,0];
edge=[1,4096];

for jj=1:nsets
    sds_id= hdfsd('select',sd_id,jj-1); % sds_id is zero-based
    fred=hdfsd('readdata',sds_id,start,[],edge);  
   spectra(:,jj)=double(fred(1:range));
end
status=hdfsd('end',sd_id);



    