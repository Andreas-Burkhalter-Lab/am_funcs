%Need to test activity vs velocity and activity vs environment

%will have 3 F-files for each plane

%cp code here to go through each

%code here for place cells

%here for df/f in f vs n
%make sure cells have the same index

%play mean F in each loc, tile 16 bars
set(0,'DefaultFigureColormap',jet)

%% Load and declare
% clear all; close all;
% function interneuron_tests(MouseID, planes, paths, name)
saveDir='G:\MA\Interneurons\';
Fs=input('Fs? ');

planes=input('Number of planes: ');
% planes=4;
MouseID = input('Mouse ID and Day: ');
% ckplanes = input('Check Plane f? (Warning slow)');
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
    %     F=zscore(F);
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
        %         forwardvel=decimate((forwardvel((1+del):(ratio*size(F,1)+del))),ratio);
        forwardvel=mean(reshape(forwardvel(1:(length(F)*ratio)),ratio,round(length(F))),1)';
        %         rotationvel=decimate((rotationvel((1+del):(ratio*size(F,1)+del))),ratio);
        rotationvel=mean(reshape(rotationvel(1:(length(F)*ratio)),ratio,round(length(F))),1)';
        
        ybinned=decimate(smooth(ybinned(1:(ratio*size(F,1))),51),ratio);
    end
    path=(path-min(path));
    %     path=path/max(path)*1e6;
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

if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
% Fall(del_cells)=[];
times=linspace(1,length(ybinned)/Fs,length(ybinned));
% Fall=Fall(:,[8,9,10,11,13,14]);
% Fall=Fall(:,[1 2 3 4 6 7 9 10 11]);
num_locs=max(pos);
num_cells=size(Fall,2);
%% Scatter dF and speed metrics
speedy=[0;diff(ybinned)];
dFcorr(num_cells,2)=0;
for cell=1:num_cells
    [ctemp,ptemp]=corrcoef((Fall(:,cell)),(forwardvel));
    dFcorr(cell,:)=[ctemp(2,1),ptemp(2,1)];
end

num_square=9;
num_figs=ceil(num_cells/num_square);
Fallt=bsxfun(@minus,Fall,min(Fall,[],1));
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    %
    %     %psd
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(speedy)),log(Fallt(:,cell)),'b+')
        title(['F vs VR speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FVR scatter', num2str(f), '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(forwardvel)),log(Fallt(:,cell)),'b+')
        title(['F vs  Forward Ball Speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FBall scatter', num2str(f), '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(forwardvel)+abs(rotationvel)),log(Fallt(:,cell)),'b+')
        title(['F vs  Total Ball Speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FTotalBall scatter', num2str(f), '.jpg']);
    
end

%% Plot correlations between dF and ball forward V

figure
hold on
corrs=dFcorr(:,1);
ps=dFcorr(:,2);

bar(corrs.*double(ps<.0005),'g')
bar(corrs.*double(ps>.0005),'r')
set(gca,'XTick',1:num_cells)
xlabel('Cell Number')
ylabel('Correlation')
title('Forward Velocity and dF/F')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'speed_dF', '.jpg']);

win=hanning(2^(floor(log2(length(Fallt)))-1));
nover=2^(floor(log2(length(Fallt)))-2);
nff=2^(floor(log2(length(Fallt)))-1);
%% Analyze speed - use diff(ybinned) and sum(forv yawv) compare with activity
figure
speedy=[0;diff(ybinned)];
speedv=forwardvel;%+rotationvel;
plot(times,speedy/max(speedy),times,speedv/max(speedv)+1);
legend('VR Speed','Ball forward+yaw')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR and Ball Speed','.jpg']);

figure
pwelch(speedy,win,nover,nff,Fs)
title('PSD Virtual reality speed')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR PSD','.jpg']);



figure
pwelch(speedv,win,nover,nff,Fs)
title('PSD Ball speed')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Ball PSD','.jpg']);



figure
pwelch(ybinned,win,nover,nff,Fs)
title('PSD Virtual reality position')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Position PSD','.jpg']);


%% xcov with speed

lags=zeros(2*length(Fall)-1,num_cells);
cross_corrs=zeros(2*length(Fall)-1,num_cells);

for cell=1:num_cells
    [cross_corrs(:,cell),lags(:,cell)]=xcov(Fall(:,cell),speedv,'coeff');
end
figure;
imagesc(cross_corrs')
title(['Cross Corr between cells and ball speed for Mouse: ', MouseID])
colormap('jet')

% figure;
% surf(cross_corrs')

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'xcov with speedv','.jpg']);

%% xcov cells
crossbox=zeros(2*length(Fall)-1,num_cells,num_cells);
lagbox=zeros(2*length(Fall)-1,num_cells,num_cells);
for cell1=1:num_cells
    for cell2=1:num_cells
        [crossbox(:,cell1,cell2),lagbox(:,cell1,cell2)]=xcorr(Fall(:,cell1),Fall(:,cell2),'coeff');
    end
end
num_square=9;
num_figs=ceil(num_cells/num_square);
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        imagesc(squeeze(crossbox(:,cell,:))');
        xlim([0 length(crossbox)])
        ylim([1 num_cells])
        title(['XCov: ', num2str(cell),' mouse ', MouseID]);
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'XCov Cells', num2str(f), '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'XCov Cells', num2str(f),'.fig']);
    
end
%% Plot each cell in 2s bins
% timelook=[1,16];
% for tl=timelook
% Fspan=round(Fs)*tl;
% vislocs=(Fspan*2):(4*Fspan):length(Fall);
% clear visPhase meanvisPhase semvisPhase
% for cell=1:num_cells
%     for r=1:(length(Fall)/(Fspan*4))
%         if vislocs(r)<=Fspan*2
%             visPhase(:,r,cell)=[zeros(Fspan*2-vislocs(r)+1,1);Fall(1:(vislocs(r)+Fspan*2),cell)];
%         elseif vislocs(r)>=(length(Fall)-Fspan*2)
%             visPhase(:,r,cell)=[Fall((vislocs(r)-Fspan*2):end,cell);zeros(vislocs(r)+Fspan*2-length(Fall),1)];
%         else
%             visPhase(:,r,cell)=Fall((vislocs(r)-Fspan*2):(vislocs(r)+Fspan*2),cell);
%         end
%     end
%     meanvisPhase(:,cell)=mean(squeeze(visPhase(:,:,cell)),2);
%     semvisPhase(:,cell)=std(squeeze(visPhase(:,:,cell)),0,2)/sqrt(size(squeeze(visPhase(:,:,cell)),2));
%     num_square=9;
% end
% num_figs=ceil(num_cells/num_square);
% for f=1:num_figs
%     fig_start=(f-1)*num_square+1;
%     if f<num_figs
%         fig_end=fig_start+(num_square-1);
%     else fig_end = num_cells;
%     end
%
%
%       figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%     splot=1;
%     for cell=fig_start:fig_end
%         subplot(sqrt(num_square),sqrt(num_square),splot)
%         hold on
%         imagesc(squeeze(visPhase(:,:,cell))')
%         xlim([1 size(visPhase,1)])
%         ylim([1 size(visPhase,2)])
%         set(gca,'XTick',[0:Fspan:size(visPhase,1)],'XTickLabel',round([-size(visPhase,1)/2:Fspan:size(visPhase,1)/2]/Fs));
%         xlabel('Seconds')
%         ylabel('Bin')
%         title(['visPhase heatmap cell: ', num2str(cell),' mouse ', MouseID]);
%         splot=splot+1;
%
%     end
%
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'visPhase heatmap', num2str(f), '.jpg']);
%     savefig([saveDir,'\',MouseID,'\', MouseID, 'visPhase heatmap', num2str(f),'.fig']);
%
%     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%     splot=1;
%     for cell=fig_start:fig_end
%         subplot(sqrt(num_square),sqrt(num_square),splot)
%         hold on
%         errorbar(meanvisPhase(:,cell),semvisPhase(:,cell))
%         line([2*Fspan+1,2*Fspan+1],[min(meanvisPhase(:,cell)-3*semvisPhase(:,cell)),max(meanvisPhase(:,cell)+2*semvisPhase(:,cell))]);
%         set(gca,'XTick',[0:Fspan:size(visPhase,1)],'XTickLabel',round([-size(visPhase,1)/2:Fspan:size(visPhase,1)/2]/Fs));
%         xlabel('Seconds')
%         xlim([1 length(meanvisPhase)])
%         title(['visPhase errorbar cell: ', num2str(cell),' mouse ', MouseID]);
%         splot=splot+1;
%
%     end
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'visPhase errorbar', num2str(f), '.jpg']);
%     savefig([saveDir,'\',MouseID,'\', MouseID, 'visPhase errorbar', num2str(f),'.fig']);
%
%
%
%
%     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%     splot=1;
%     %plot all F with lines to show split
%     for cell=fig_start:fig_end
%         subplot(sqrt(num_square),sqrt(num_square),splot)
%         t=0:(size(meanvisPhase,1)-1);
%         plot(t,meanvisPhase(:,cell),'b',t,(meanvisPhase(:,cell)+2*semvisPhase(:,cell)),'g:',t,(meanvisPhase(:,cell)-2*semvisPhase(:,cell)),'g:')
%         set(gca,'XTick',[t(1:Fspan:end)],'XTickLabel',round([-size(meanvisPhase,1)/2:Fspan:size(meanvisPhase,1)/2]/Fs));
%         xlabel('Seconds')
%         xlim([0, length(t)])
%         line([2*Fspan,2*Fspan],[min(meanvisPhase(:,cell)-3*semvisPhase(:,cell)),max(meanvisPhase(:,cell)+2*semvisPhase(:,cell))]);
%         title(['visPhase confidence cell: ', num2str(cell),' mouse ',MouseID]);
%         splot=splot+1;
%
%     end
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'visPhase confidence', num2str(f), '.jpg']);
%     savefig([saveDir,'\',MouseID,'\', MouseID, 'visPhase confidence', num2str(f),'.fig']);
%
%
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         [~,prtainds]=max(meanvisPhase);
%         [~,sortedinds]=sort(prtainds);
%         imagesc(meanvisPhase(:,sortedinds)')
%         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'visPhase heatmap all', num2str(f), '.jpg']);
%         savefig([saveDir,'\',MouseID,'\', MouseID, 'visPhase heatmap all', num2str(f),'.fig']);
%
%
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         surf(meanvisPhase(:,sortedinds)')
%
% end
%
% end
%


%% Look at PRTA

[PRTA,meanPRTA,semPRTA]=calc_PRTA(Fall,rewards,Fs,MouseID,saveDir);


%% power spectrum
[pwn,fpwn]=pwelch(Fall(:,:),win,nover,nff,Fs);

%% plotting
num_square=9;
num_figs=ceil(num_cells/num_square);
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    
    %psd
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    suptitle(['PSD for cells: ',num2str(fig_start),' to ',num2str(fig_end), ' Mouse: ',MouseID])
    
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        len=round(length(fpwn)/10);
        plot(fpwn,(10*log10(pwn(:,cell))))
        xlim([0 .05])
        if ismember(cell,npil)
            title('Neurites')
        elseif ismember(cell, nothing)
            title('Nothing Visible')
        elseif ismember(cell,missing)
            title(['(Missing) Cell ', num2str(cell)])
        else
            title(['Cell ', num2str(cell)])
        end
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD_small', num2str(f), '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD_small', num2str(f),'.fig']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    suptitle(['PSD for cells: ',num2str(fig_start),' to ',num2str(fig_end), ' Mouse: ',MouseID])
    
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        len=round(length(fpwn)/10);
        plot(fpwn,(10*log10(pwn(:,cell))))
        if ismember(cell,npil)
            title('Neurites')
        elseif ismember(cell, nothing)
            title('Nothing Visible')
        elseif ismember(cell,missing)
            title(['(Missing) Cell ', num2str(cell)])
        else
            title(['Cell ', num2str(cell)])
        end
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f), '.svg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f),'.fig']);
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    suptitle(['F traces for cells: ',num2str(fig_start),' to ',num2str(fig_end), ' Mouse: ',MouseID])
    splot=1;
    %plot all F with lines to show split
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        cellF=Fall(:,cell);
        plot(cellF)
        xlim([0 length(Fall)])
        if ismember(cell,npil)
            title('Neurites')
        elseif ismember(cell, nothing)
            title('Nothing Visible')
        elseif ismember(cell,missing)
            title(['(Missing) Cell ', num2str(cell)])
        else
            title(['Cell ', num2str(cell)])
        end
        splot=splot+1;
        set(gca,'Fontsize',18)
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenv', num2str(f), '.svg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenv', num2str(f),'.fig']);
    
end

%% connection test
order=1:num_cells;
% order=[1 3 21 2 8 18 19 22 4 7 9 17 20 12 13];
% sorder=strsplit(num2str([order, npil, nothing]));
figure
% imagesc(corrcoef(Fall(:,[order, npil, nothing])))
title(['Correlation Matrix for mouse: ',MouseID])
imagesc(corrcoef(Fall))

colorbar
% set(gca,'Xtick',1:length(sorder),'XTickLabel',sorder)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'CorrMat', '.jpg']);

%% Filter freq bands
lsig=zeros(size(Fall));
hsig=lsig;
gsig=hsig;
asig=gsig;
betsig=asig;
thesig=betsig;
delsig=thesig;
for cell=1:num_cells
    lp=designfilt('lowpassfir','FilterOrder',32, ...
        'PassbandFrequency',.2,'StopbandFrequency',.4,...
        'SampleRate',Fs);
    lsig(:,cell)=filtfilt(lp,Fall(:,cell));
    del=mean(grpdelay(lp));
    hp=designfilt('highpassfir','FilterOrder',16, ...
        'PassbandFrequency',.2,'StopbandFrequency',.1,...
        'SampleRate',Fs);
    if Fs>=6
        delp=designfilt('bandpassfir','FilterOrder',16, ...
            'CutoffFrequency1',.2,'CutoffFrequency2',3,...
            'SampleRate',Fs);
        delsig(:,cell)=filtfilt(delp,Fall(:,cell));
    end
    if Fs>=16
        thep=designfilt('bandpassfir','FilterOrder',16, ...
            'CutoffFrequency1',4,'CutoffFrequency2',8,...
            'SampleRate',Fs);
        thesig(:,cell)=filtfilt(thep,Fall(:,cell));
    end
    if Fs>=30
        ap=designfilt('bandpassfir','FilterOrder',16, ...
            'CutoffFrequency1',8,'CutoffFrequency2',15,...
            'SampleRate',Fs);
        asig(:,cell)=filtfilt(ap,Fall(:,cell));
    end
    if Fs>=60
        betp=designfilt('bandpassfir','FilterOrder',16, ...
            'CutoffFrequency1',16,'CutoffFrequency2',30,...
            'SampleRate',Fs);
        betsig(:,cell)=filtfilt(betp,Fall(:,cell));
    end
    if Fs>=75
        gp=designfilt('highpassfir','FilterOrder',16, ...
            'PassbandFrequency',30,'StopbandFrequency',25,...
            'SampleRate',Fs);
        gsig(:,cell)=filtfilt(gp,Fall(:,cell));
    end
    
    hsig(:,cell)=filtfilt(hp,Fall(:,cell));
end
clear angle
% lsig, gsig, betsig,thesig,delsig,hsig
[PRTA,meanPRTA,semPRTA]=calc_PRTA(lsig,rewards,Fs,[MouseID,'lowfreq'],saveDir);
[PRTA,meanPRTA,semPRTA]=calc_PRTA(angle(hilbert(lsig))/pi,rewards,Fs,[MouseID,'lowfreq phase'],saveDir);

[PRTA,meanPRTA,semPRTA]=calc_PRTA(delsig,rewards,Fs,[MouseID,'deltafreq'],saveDir);
[PRTA,meanPRTA,semPRTA]=calc_PRTA(angle(hilbert(delsig)),rewards,Fs,[MouseID,'deltafreq phase'],saveDir);

[PRTA,meanPRTA,semPRTA]=calc_PRTA(abs(hilbert(thesig)),rewards,Fs,[MouseID,'thetafreq'],saveDir);
[PRTA,meanPRTA,semPRTA]=calc_PRTA(angle(hilbert(thesig)),rewards,Fs,[MouseID,'thetafreq phase'],saveDir);

[PRTA,meanPRTA,semPRTA]=calc_PRTA(abs(hilbert(betsig)),rewards,Fs,[MouseID,'betafreq'],saveDir);
[PRTA,meanPRTA,semPRTA]=calc_PRTA(angle(hilbert(betsig)),rewards,Fs,[MouseID,'betafreq phase'],saveDir);

[PRTA,meanPRTA,semPRTA]=calc_PRTA(abs(hilbert(asig)),rewards,Fs,[MouseID,'alphafreq'],saveDir);
[PRTA,meanPRTA,semPRTA]=calc_PRTA(angle(hilbert(asig)),rewards,Fs,[MouseID,'alphafreq phase'],saveDir);

[PRTA,meanPRTA,semPRTA]=calc_PRTA(abs(hilbert(gsig)),rewards,Fs,[MouseID,'gammafreq'],saveDir);
[PRTA,meanPRTA,semPRTA]=calc_PRTA(angle(hilbert(gsig)),rewards,Fs,[MouseID,'gammafreq phase'],saveDir);
close all

%% More freq analysis
params.fpass = [0 Fs/2];
params.Fs=Fs;
params.tapers=[5 9];
movingwin = [5 1];
for jj=1:num_cells
    % jj=1;
    %Spectrogram
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    suptitle(['Spectrogram for cell: ',num2str(jj), ' Mouse: ',MouseID])
    
    subplot(4,1,1:2)
    [specgram_data.s,specgram_data.t,specgram_data.f]=mtspecgramc(Fall(:,jj),movingwin,params);
    imagesc(specgram_data.t,specgram_data.f, log(specgram_data.s'))
    
    
    subplot(4,1,3:4)
    hold on
    plot(times,pos/max(pos)+(max(Fall(:,jj))))
    
    plot(times,speedv/max(speedv)+(min(Fall(:,jj))-1))
    plot(times,Fall(:,jj));
    xlim([0 max(times)])
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'spectrogram Tres', num2str(jj), '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Tres', num2str(jj),'.fig']);
end
movingwin = [50 20];
for jj=1:num_cells
    % jj=1;
    %Spectrogram
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    suptitle(['Spectrogram for cell: ',num2str(jj), ' Mouse: ',MouseID])
    
    subplot(4,1,1:2)
    [specgram_data.s,specgram_data.t,specgram_data.f]=mtspecgramc(Fall(:,jj),movingwin,params);
    imagesc(specgram_data.t,specgram_data.f, log(specgram_data.s'))
    
    
    subplot(4,1,3:4)
    hold on
    plot(times,pos/max(pos)+(max(Fall(:,jj))))
    
    plot(times,speedv/max(speedv)+(min(Fall(:,jj))-1))
    plot(times,Fall(:,jj));
    xlim([0 max(times)])
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'spectrogram Fres', num2str(jj), '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Fres', num2str(jj),'.fig']);
end

%% Slow oscillations with speed, ybinned
swin=300;
ewin=420;
% swin=60;
% ewin=120;
rewsmall2=zeros(length(Fall),1);
num_cells=size(Fall,2);
for jj=1:length(Fall)
    rewsmall2(jj)=max(rewards(((jj-1)*ratio+1):(ratio*jj)));
end
[~,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs)*2);
avgrewdist=median(diff(rewlocs));
sframe=round(swin*Fs);
eframe=round(ewin*Fs);
constrain = @(sig) (sig-min(sig(sframe:eframe)))/max(sig(sframe:eframe)-min(sig(sframe:eframe)));
for cell=1:num_cells
    gmag=abs(hilbert(gsig(:,cell)));
    amag=abs(hilbert(asig(:,cell)));
    
    %         for n=2:2
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    %     offset=max(Fall(:,cell))/max(Fall(sframe:eframe,cell));
    offset=1;
    suptitle(['Oscillations for cell: ',num2str(cell), ' Mouse: ',MouseID])
    subplot(3,2,1)
    hold on
    %     plot(times,Fall(:,cell)/max(Fall(sframe:eframe,cell))+offset*2,...
    %         times,hsig(:,cell)/max(hsig(sframe:eframe,cell))+offset,...
    %         times,lsig0(:,cell)/max(lsig0(sframe:eframe,cell)));
    %     if Fs==100
    %         plot(times,gsig(:,cell)/max(gsig(sframe:eframe,cell))+3*offset);
    %     end
    %     plot(times,pos/max(pos)+min(Fall(sframe:eframe,1)-1));%,times,(Fall(:,1)))
    plot(times,constrain(Fall(:,cell))+offset*2,...
        times,constrain(hsig(:,cell))+offset,...
        times,constrain(lsig(:,cell)));
    legend('dF/F','Above .05','Below .05')
    if Fs==100
        %         plot(times,constrain(gsig(:,cell))+3*offset);
        plot(times,constrain(gmag)+3*offset);
        plot(times,constrain(amag)+4*offset);
        legend('dF/F','Above .05','Below .05','Gamma','Beta')
    end
    plot(times,pos/max(pos)-1);%,times,(Fall(:,1)))
    xlim([swin ewin])
    
    title('Activity and Position')
    subplot(3,2,3)
    hold on
    plot(times,constrain(Fall(:,cell))+offset*2,...
        times,constrain(hsig(:,cell))+offset,...
        times,constrain(lsig(:,cell)));
    legend('dF/F','Above .05','Below .05')
    
    if Fs==100
        plot(times,constrain(gmag)+3*offset);
        plot(times,constrain(amag)+4*offset);
        legend('dF/F','Above .05','Below .05','Gamma','Beta')
        
    end
    plot(times,speedv/max(speedv)-1);%,times,(Fall(:,1)))
    title('Activity and VR speed')
    xlim([swin ewin])
    
    subplot(3,2,5)
    hold on
    plot(times,constrain(Fall(:,cell))+offset*2,...
        times,constrain(hsig(:,cell))+offset,...
        times,constrain(lsig(:,cell)));
    legend('dF/F','Above .05','Below .05')
    
    if Fs==100
        plot(times,constrain(gmag)+3*offset);
        plot(times,constrain(amag)+4*offset);
        legend('dF/F','Above .05','Below .05','Gamma','Beta')
    end
    plot(times,speedy/max(speedy)-1);%,times,(Fall(:,1)))
    title('Activity and Ball Speed')
    xlim([swin ewin])
    
    subplot(3,2,2)
    plot(fpwn,(10*log10(pwn(:,cell))))
    
    subplot(3,2,4)
    plot(fpwn,(10*log10(pwn(:,cell))))
    xlim([0 .05])
    
    
    subplot(3,2,6)
    plot(cross_corrs(:,cell))
    title('Cross Corr between cell and ball speed')
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'close look ', num2str(cell),',', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'close look', num2str(cell),',','.fig']);
    
    %Plot just df/F with speed and pos and rew
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    plot(times,(constrain(Fall(:,cell))+1))
    hold on
    plot(times,speedv/max(speedv))
    plot(times,abs(speedy)/max(speedy)-1)
    plot(times,pos/max(pos)-1)
    plot(times,rotationvel/max(rotationvel))
    plot(times,rewsmall2+1)
    xlim([swin ewin])
    legend('dF','Ball Forward','VR Speed', 'VR Pos','Rewards')
    title(['Speed and dF for Cell ', num2str(cell)])
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Pos and Speed look ', num2str(cell),',', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Pos and Speed look', num2str(cell),',','.fig']);
    
    
    %Plot dF VR pos and VR speed on one axis
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    plot(times,((Fall(:,cell))))
    hold on
    %     plot(times,speedv/max(speedv))
    plot(times,abs(speedy)/max(speedy))
    plot(times,pos/max(pos))
    %     plot(times,rotationvel/max(rotationvel))
    %     plot(times,rewsmall+1)
    xlim([swin ewin])
    ylim([min(Fall(:,cell)), max(Fall(:,cell))])
    legend('dF','VR Speed', 'VR Pos')
    title(['Speed and dF for Cell ', num2str(cell)])
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR speed same axis look ', num2str(cell),',','.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'VR speed same axis look', num2str(cell),',','.fig']);
    
    
    %         end
end



%%

%all traces overlayed
figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
for cell=1:num_cells
    
    plot(fpwn,(10*log10(pwn(:,cell))))
end
title('PSD Overlay')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
set(gca,'Fontsize',24,'Xtick',linspace(0,2,5))
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD overlay','.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD overlay','.fig']);

figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
plot(fpwn,(10*log10(mean(pwn,2))))
title('Mean PSD')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
set(gca,'Fontsize',24,'Xtick',linspace(0,2,5))

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD multi ', '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD multi', '.fig']);

figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
for cell=1:num_cells
    plot(times,Fall(:,cell))
end
xlabel('Seconds')
ylabel('dF/F')
set(gca,'Fontsize',24)

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dF overlay ', '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dF overlay', '.fig']);
figure('units','normalized', 'Position', [.01 .05 .98 .87]); hold on
for cell=1:num_cells
    plot(times,Fall(:,cell))
end
xlim([0 120])
xlabel('Seconds')
ylabel('dF/F')
set(gca,'Fontsize',24)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dF overlay 2', '.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dF overlay 2', '.fig']);
