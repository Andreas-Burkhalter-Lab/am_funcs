%% Load and declare
clear all; close all;
% function interneuron_tests(MouseID, planes, paths, name)
saveDir='F:\MA Data\Interneurons\';
Fs=input('Fs? ');

planes=input('Number of planes: ');
% planes=4;
MouseID = input('Mouse ID and Day: ');

usecell=input('Use cell: ?');
% MouseID='E31D1';
names{planes}=0;
% names={'1119_000_000_hmm1_roibyhand_F_within.mat','1119_000_000_hmm2_roibyhand_F_within.mat',...
%     '1119_000_000_hmm3_roibyhand_F_within.mat','1119_000_000_hmm4_roibyhand_F_within.mat'};
paths{planes}=0;
% paths={'E:\2015\E31\151119\','E:\2015\E31\151119\',...
%     'E:\2015\E31\151119\','E:\2015\E31\151119\'};
% num_files=1;

npil=[];
nothing=[];
missing=[];
% nothing=[15 20 25]; %E46
% nothing=[4 9 10 15 22]; %X14
% missing=[9]; %E46 days 4 5
% missing=[11]; %x14 days 2 3 4
% % Fcell{planes}=0;
Fall=[];
% Features{planes}=0;
% Probs{planes}=0;

% v=[];
%%

for p=1:planes
    disp(['Plane ', num2str(p)]);
    [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    %     files{p}=name{p}(1:end-4);
    load([paths{p},names{p}]);
    lp=designfilt('lowpassfir','FilterOrder',128, ...
        'PassbandFrequency',5,'StopbandFrequency',7,...
        'SampleRate',length(forwardvel)/length(F)*Fs);
    F=Fc2;
    %     F=bsxfun(@rdivide,F,max(F,[],1));
    %     F=bsxfun(@minus,F,mean(F,1));
    %Downsample ybinned to match F
    if length(ybinned) ~= size(F,1)
        rewratio=length(ybinned)/size(F,1);
        ratio=floor(length(ybinned)/size(F,1));
        path=downsample(ybinned(1:(ratio*size(F,1))),ratio);
        
        for jj=1:200:(length(forwardvel)-200)
            forwardvel(jj:(jj+199))=max(forwardvel(jj:(jj+199)));
            rotationvel(jj:(jj+199))=max(rotationvel(jj:(jj+199)));
        end
        forwardvel=filter(lp,forwardvel);
        rotationvel=filter(lp,rotationvel);
        del=mean(grpdelay(lp));
        %         (1+del):(ratio*size(F,1)+del))49:ratio*size(F,1)+48
        %         forwardvel=decimate((forwardvel((1+del):(ratio*size(F,1)+del)))),ratio);
        forwardvel=mean(reshape(forwardvel(1:(length(F)*ratio)),ratio,round(length(F))),1)';
        %         rotationvel=decimate((rotationvel((1+del):(ratio*size(F,1)+del))),ratio);
        rotationvel=mean(reshape(rotationvel(1:(length(F)*ratio)),ratio,round(length(F))),1)';
        ybinned=decimate(smooth(ybinned(1:(ratio*size(F,1))),51),ratio);
    end
    path=(path-min(path));
    path=path/max(path)*180/5;
    pos=ceil(path+eps);
    if ~isempty(Fall)
        if length(Fall)<length(F)
            Fall=[Fall,F(1:length(Fall),:)];
        else
            Fall=[Fall(1:length(F),:), F];
        end
    else
        Fall=F;
    end
    %     novel_start=round(timeSplit(1)/rewratio);
    %     novel_end=round(timeSplit(2)/rewratio);
end

if ~exist([saveDir,'\',MouseID,'\'],'dir')
    mkdir([saveDir,'\',MouseID,'\']);
end
% Fall(del_cells)=[];
times=linspace(1,length(ybinned)/Fs,length(ybinned));

num_locs=max(pos);
num_cells=size(Fall,2);

win=hanning(2^(floor(log2(length(Fall)))-1));
nover=2^(floor(log2(length(Fall)))-2);
nff=2^(floor(log2(length(Fall)))-1);
% usecell=9;
celltrace=Fall(:,usecell);
%% figure one - dF and ball speed cell 1
speedy=[0;diff(ybinned)];
speedv=forwardvel;%+rotationvel;
speedsum=forwardvel+rotationvel;
period=[0,120;320,440;960,1080];

for s=1:3
figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
plot(times,celltrace,'linewidth',1.5);
plot(times,speedsum/max(speedsum)+min(celltrace)-1,'linewidth',1.5)
line([0,0],[-.4,.6],'Color','k','linewidth',3)
legend('dF/F','Running Speed')
xlim(period(s,:))
xlabel('Seconds')
% ylabel('Running Speed                                                            1 dF/F                                 ')
set(gca,'YTick',[],'fontsize',24)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dF and speed cell',num2str(usecell),' period ' num2str(s),'.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dF and speed cell',num2str(usecell),' period ' num2str(s),'.fig']);
end

%% Figure Two

speedy=[0;diff(ybinned)];
speedv=forwardvel;%+rotationvel;
speedsum=forwardvel+rotationvel;

figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
plot(times,celltrace,'linewidth',1.5,'Color','b');
plot(times,speedsum/max(speedsum)+min(celltrace)-1,'linewidth',1.5,'Color','r')

plot(times,celltrace-max(celltrace)-1,'linewidth',1.5,'Color','b');
plot(times,(speedsum/max(speedsum))...
    +min(celltrace)-1-max(celltrace),'linewidth',1.5,'color','r')
legend('dF/F','Running Speed')
line([0,0],[-.4,.6],'Color','k','linewidth',3)
line([0,0],[-2.8,-1.8],'Color','k','linewidth',3)

xlim([0,120])
xlabel('Seconds')
% ylabel('1 dF/F                                      Running Speed                           1 dF/F                 ')
set(gca,'YTick',[],'fontsize',24)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dF and speed Overlay cell ',num2str(usecell), '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dF and speed Overlay cell ',num2str(usecell), '.fig']);



%% Figure 3 scatter, lets look at linear/log VR/ball

figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
s1=subplot(2,1,1);

scatter((abs(speedy))/max(abs(speedy)),celltrace,'b.')
lsline(s1)
r=corrcoef(abs(speedy),celltrace);
text(.9,1,['r=',num2str(r(2,1))])
title(['F vs VR speed'])
ylabel('dF/F')
xlabel('VR Speed')
set(gca,'fontsize',24,'Xtick',linspace(0,1,5))

% s2=subplot(2,2,2);
% scatter(log(abs(speedy)),log(celltrace-min(celltrace)),'b.')
% h=lsline(s2);
% set(h(1),'color','r')
% r=corrcoef(log(abs(speedy)),log((celltrace-min(celltrace))));
% text(.9,1,['r=',num2str(r(2,1))])
% title(['F vs VR speed, Log Scale'])
% ylabel('dF/F')
% xlabel('VR Speed')


s3=subplot(2,1,2);

scatter((abs(speedsum))/max(abs(speedsum)),celltrace,'b.')
lsline(s3)
r=corrcoef(abs(speedsum),celltrace);
text(.9,1,['r=',num2str(r(2,1))])
title(['F vs Running Speed'])
ylabel('dF/F')
xlabel('Running Speed')
set(gca,'fontsize',24,'Xtick',linspace(0,1,5))
% s4=subplot(2,2,4);
% scatter(log(abs(speedsum)),log(celltrace-min(celltrace)),'b.')
% h=lsline(s4);
% set(h(1),'color','r')
% r=corrcoef(log(abs(speedsum)),log((celltrace-min(celltrace))));
% text(.9,1,['r=',num2str(r(2,1))])
% title(['F vs Running speed, Log Scale'])
% ylabel('dF/F')
% xlabel('Running Speed')

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'scatter cell ',num2str(usecell), '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'scatter cell ',num2str(usecell), '.fig']);




%% PSD

[pwn,fpwn]=pwelch(celltrace,win,nover,nff,Fs);


figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
subplot(2,1,1)
plot(fpwn,(10*log10(pwn)))
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
set(gca,'fontsize',24,'Xtick',linspace(0,2,5))

subplot(2,1,2)
plot((fpwn),(10*log10(pwn)))
xlim([0 .5])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
set(gca,'fontsize',24,'Xtick',linspace(0,.5,5))

% suptitle('PSD of dF/F')

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD cell ',num2str(usecell), '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD cell ',num2str(usecell), '.fig']);




%% add rewards indicators


rewsmall2=zeros(size(speedy));
for jj=1:length(speedy)
    rewsmall2(jj)=max(rewards(((jj-1)*ratio+1):(ratio*jj)));
end
[rewsmall,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs)*2);
avgrewdist=median(diff(rewlocs));

figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
plot(times,celltrace,'linewidth',1.5);
plot(times,speedsum/max(speedsum)+min(celltrace)-1,'linewidth',1.5)
for r=rewlocs
    line([rewlocs/Fs rewlocs/Fs],[min(celltrace)-1 max(celltrace)],'color','k')
end
legend('dF/F','Running Speed','Rewards')
line([0,0],[-.4,.6],'Color','k','linewidth',3)
xlim([0,120])
xlabel('Seconds')
ylabel(' Running Speed                                                       1 dF/F                                                                   ')
set(gca,'YTick',[])
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dF and speed rewards cell ',num2str(usecell), '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dF and speed rewards cell ',num2str(usecell), '.fig']);

%% plot with VR pos not ball speed

figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
plot(times,celltrace,'linewidth',1.5);
plot(times,pos/max(pos)+min(celltrace)-1,'linewidth',1.5)
for r=rewlocs
    line([rewlocs/Fs rewlocs/Fs],[min(celltrace)-1 max(celltrace)],'color','k')
end
legend('dF/F','VR position','Rewards')
line([0,0],[-.4,.6],'Color','k','linewidth',3)
xlim([0,120])
xlabel('Seconds')
% ylabel('VR Position                                                         1 dF/F                                    ')
set(gca,'YTick',[],'fontsize',24)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dF and pos rewards cell ',num2str(usecell), '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dF and pos rewards cell ',num2str(usecell), '.fig']);

%% extend time



figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
subplot(2,1,1)
plot(times,celltrace,'linewidth',1.5);
for r=rewlocs
    line([rewlocs/Fs rewlocs/Fs],[min(celltrace) max(celltrace)],'color','k')
end
legend('dF/F','Rewards')
line([0,0],[-.4,.6],'Color','k','linewidth',3)
xlim([0,120])
xlabel('Seconds')
% ylabel('1 dF/F                      ')
set(gca,'YTick',[],'fontsize',24)
subplot(2,1,2)
plot(times,celltrace,'linewidth',1.5);
for r=rewlocs
    line([rewlocs/Fs rewlocs/Fs],[min(celltrace) max(celltrace)],'color','k')
end
line([120,120],[-.4,.6],'Color','k','linewidth',3)
% ylabel('1 dF/F                      ')
set(gca,'YTick',[],'fontsize',24)

xlim([120,240])
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dF extended rewards cell ',num2str(usecell), '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dF extended rewards cell ',num2str(usecell), '.fig']);

%% PRTA
timelook=[1, round(avgrewdist/Fs*2/3), 16];
for tl=timelook
    Fspan=round(Fs)*tl;
    % rewlocs=find(rews);
    PRTA = zeros(Fspan*4+1,length(rewlocs));
    meanPRTA=zeros(1,Fspan*4+1);
    semPRTA=meanPRTA;
    for r=1:length(rewlocs)
        if rewlocs(r)<=Fspan*2
            PRTA(:,r)=[zeros(Fspan*2-rewlocs(r)+1,1);celltrace(1:(rewlocs(r)+Fspan*2))];
        elseif rewlocs(r)>=(length(celltrace)-Fspan*2)
            %             PRTA(:,r,cell)=[Fall((rewlocs(r)-Fspan*2):end,cell);zeros(rewlocs(r)+Fspan*2-length(rewsmall),1)];
            PRTA(:,r)=[celltrace((rewlocs(r)-Fspan*2):end);zeros(size(PRTA,1)-length(celltrace((rewlocs(r)-Fspan*2):end)),1)];
        else
            PRTA(:,r)=celltrace((rewlocs(r)-Fspan*2):(rewlocs(r)+Fspan*2));
        end
    end
    meanPRTA(:)=mean(squeeze(PRTA(:,:)),2);
    semPRTA(:)=std(squeeze(PRTA(:,:)),0,2)/sqrt(size(squeeze(PRTA(:,:)),2));
    
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    imagesc(PRTA')
    set(gca,'Ydir','normal')
    xlim([1 size(PRTA,1)])
    set(gca,'XTick',[0:(Fs*tl):size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs),'fontsize',24)
    xlabel('Seconds')
    ylim([1 size(PRTA,2)])
    line([2*Fspan,2*Fspan], [1 size(PRTA,2)],'linewidth',2,'color','r')
    title(['Reward Activity heatmap']);
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap ',num2str(tl*4),' cell ',num2str(usecell), '.svg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap ',num2str(tl*4),' cell ',num2str(usecell), '.fig']);
    
    
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
    t=0:(size(meanPRTA,2)-1);
    plot(t,meanPRTA(:),'b',t,(meanPRTA+2*semPRTA),'r--',t,(meanPRTA-2*semPRTA),'r--')
    set(gca,'XTick',[t(1:Fspan:end)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs),'fontsize',24);
    xlabel('Seconds')
    line([2*Fspan,2*Fspan],[min(meanPRTA(:)-3*semPRTA(:)),max(meanPRTA(:)+2*semPRTA(:))]);
    line([2*Fspan-avgrewdist,2*Fspan-avgrewdist],[min(meanPRTA(:)-3*semPRTA(:)),max(meanPRTA(:)+2*semPRTA(:))],'Color','r');
    line([2*Fspan+avgrewdist,2*Fspan+avgrewdist],[min(meanPRTA(:)-3*semPRTA(:)),max(meanPRTA(:)+2*semPRTA(:))],'Color','r');
    xlim([0 size(PRTA,1)])
    title(['Mean Reward Averages']);
    
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA confidence cell ',num2str(usecell),' '  num2str(tl*4),' seconds', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA confidence cell ',num2str(usecell),' ', num2str(tl*4),' seconds','.fig']);
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
    subplot(2,2,1)
    plot(fpwn,(10*log10(pwn)))
    set(gca,'YTick',[],'fontsize',24)
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    
    
    subplot(2,2,2)
    imagesc(PRTA')
    %     set(gca,,'Ydir','normal')
        line([2*Fspan,2*Fspan], [1 size(PRTA,2)],'linewidth',2,'color','r')
    xlim([1 size(PRTA,1)])
    set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs),'fontsize',24,'Ydir','normal');
    xlabel('Seconds')
    ylim([1 size(PRTA,2)])
    title(['Reward Activity heatmap']);
  
    subplot(2,2,3:4)
    plot(times,celltrace,'linewidth',1.5);
    for r=rewlocs
        line([rewlocs/Fs rewlocs/Fs],[min(celltrace) max(celltrace)],'color','k')
    end
    legend('dF/F','Rewards')
    line([0,0],[-.4,.6],'Color','k','linewidth',3)
    xlim([0,240])
    text(210,1,['Avg Reward Frequency=',num2str(1/(avgrewdist/Fs))])
    xlabel('Seconds')
    ylabel('1 dF/F                      ')
    set(gca,'YTick',[],'fontsize',24)
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD and PRTA heatmap cell ',num2str(usecell),' '  num2str(tl*4),' seconds', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD and PRTA heatmap cell ',num2str(usecell),' ', num2str(tl*4),' seconds','.fig']);
    
  
end

%% Filter at .3Hz

lp = designfilt('lowpassfir', 'PassbandFrequency', .3, 'StopbandFrequency',...
    .55, 'PassbandRipple', .001, 'StopbandAttenuation', 60, 'SampleRate', 15.5/4);
lsig=filter(lp,celltrace);
ldel=mean(grpdelay(lp));
figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
hold on
plot(times,celltrace+max(lsig))
plot(times(1:(end-ldel)),lsig((1+ldel):end))
plot(times,pos/max(pos)-1+min(lsig));
% for r=rewlocs
%     line([rewlocs/Fs rewlocs/Fs],[min(celltrace)-1 max(celltrace)+1.5],'color','k')
% end
legend('dF/F','Lowpass at .5Hz','Running Speed')
line([0,0],[-.4,.6],'Color','k','linewidth',3)
xlim([0,120])
xlabel('Seconds')
% ylabel(' Running Speed                                                       1 dF/F                                                    ')
set(gca,'YTick',[],'fontsize',24)

   saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'lowpass trace cell ',num2str(usecell), '.svg']);
  savefig([saveDir,'\',MouseID,'\', MouseID, 'lowpass trace cell ',num2str(usecell),'.fig']);
 












