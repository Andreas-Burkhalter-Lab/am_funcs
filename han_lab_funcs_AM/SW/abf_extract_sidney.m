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



