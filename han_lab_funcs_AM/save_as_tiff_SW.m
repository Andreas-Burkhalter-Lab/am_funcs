
[fname,pname]=uigetfile('.mat','Choose Video File');
first=input('First Frame: ');
last=input('Last Frame: ');
imgname=input('Name Video: ');
imgname=[imgname,'.tif'];
load([pname,fname]);
for k=1:(last-first)
    imwrite(uint16(squeeze(video(:,:,(k+first-1)))),imgname,'WriteMode','append','Compression','none');
end
