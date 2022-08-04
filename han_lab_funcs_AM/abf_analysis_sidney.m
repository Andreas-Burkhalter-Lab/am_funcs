clear all;
%Choose abf file to load
[filename,path]=uigetfile('*.abf','pick your file');
cd (path); %set path
fullfilename=[path char(filename)];
data=abfload(fullfilename);

timebin=5*60*1000;
datamat=zeros(timebin,(floor(length(data)/(timebin))),size(data,2));
for d=1:8
    for i=1:6
        temp=data(((i-1)*timebin+1):(i*timebin),d);
        datamat(:,i,d)=temp;
    end
end
figure('units','normalized', 'Position', [.01 .05 .65 .43]);
for f=1:6
    subplot(3,2,f)
    plot(squeeze(datamat(:,f,2)))
    title(['Position in Time Bin', num2str(f)])
end
figure('units','normalized', 'Position', [.01 .05 .65 .43]);
for f=1:6
    subplot(3,2,f)
    forward=squeeze(datamat(:,f,6));
    rotation=squeeze(datamat(:,f,7));
    plot(sqrt(forward.^2+rotation.^2))
    title(['Total Velocity in Time Bin', num2str(f)])
end 

%Look at locomotor acitivity between sessions
datacell{2}=1;
total{2}=1;
figure('units','normalized', 'Position', [.01 .05 .65 .43]);
for ii=1:2
[filename,path]=uigetfile('*.abf','pick your file');
cd (path); %set path
fullfilename=[path char(filename)];
timebin=5*60*1000;
datacell{ii}=abfload(fullfilename);
total{ii}=sqrt(datacell{ii}(:,6).^2+datacell{ii}(:,7).^2);
histogram(total{ii})
hold on

end
%Change Legend Entries for histograms here
legend({'First Session','Second Session'})
%Statistical tests here
[httest,pttest]=ttest2(total{1},total{2})

[pranksum,hranksum,stats]=ranksum(total{1},total{2},'tail','right')


%     angle=data(scanstart:scanstop,4); %anglecut=view angle, channel 4 of "data"
% rewards=data(scanstart:scanstop,1);  %rewcut=rewards
% galvo=data(scanstart:scanstop,5);    %131018
% forwardvel=data(scanstart:scanstop,6);   %120806 EH
% rotationvel=data(scanstart:scanstop,7);  %120806 EH
% % plot(rewcut);ginput(1);
% ybinned=data(scanstart:scanstop,2);    %Ycut=yposition
% numframes=length(scanstart:scanstop);
% if size(data,2)>7%131018
%     ch8=data(scanstart:scanstop,8); %depending on pclamp protocol, may or may not have ch8
% end