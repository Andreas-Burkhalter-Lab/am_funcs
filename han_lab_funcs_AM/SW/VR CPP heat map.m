clear all;
[filename,path]=uigetfile('*.abf','pick your file');
cd (path); %set path
fullfilename=[path char(filename)];
data=abfload(fullfilename);

angle=data(:,4); %anglecut=view angle, channel 4 of "data"
rewards=data(:,1);  %rewcut=rewards
galvo=data(:,5);    %131018
forwardvel=data(:,6);   %120806 EH
rotationvel=data(:,7);  %120806 EH
% plot(rewcut);ginput(1);
ybinned=data(:,2);    %Ycut=yposition
xbinned=data(:,3);


ybinned=ybinned-min(ybinned(:));
ybinned=ybinned/max(ybinned)*180;
figure;
ybinnedtemp=floor(ybinned);
tempinds=logical(ybinned<179) & logical(ybinned>1);
histogram(ybinned,180);
figure;
subplot(8,1,1:7)
numbins=180;
xbins=numbins/2;
n=histcounts(ybinned(tempinds),numbins);
tempmat=repmat(n',1,xbins);
cmap=colormap('parula');

imshow(tempmat/max(tempmat(:))*65,cmap);
set(gca,'XTick',1:50)
% colorbar
hold on
xbint=xbinned-min(xbinned);
xbint=xbint/max(xbint)*xbins;
ybint=ybinned/max(ybinned);
ybint=ybint*numbins;
plot(xbint,ybint,'r')

subplot(8,1,8)
histogram(ybinned(tempinds),numbins);