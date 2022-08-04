%Need to test activity vs velocity and activity vs environment

%will have 3 F-files for each plane

%cp code here to go through each

%code here for place cells

%here for df/f in f vs n
%make sure cells have the same index

%play mean F in each loc, tile 16 bars

%% Load and declare
% clear all; close all;
function [Fall,rotation,forward]=TwoEnvInterneuronTest(varargin)
set(0,'DefaultFigureColormap',jet)

saveDir='F:\MA Data\Interneurons\';
if nargin<1
    Fs=input('Fs? ');
    planes=input('Number of planes: ');
    MouseID = input('Mouse ID and Day: ');
    manual=1;
else
    Fs=varargin{1};
    planes=varargin{2};
    MouseID=varargin{3};
    manual=0;
end
paths{planes}=0;
names{planes}=0;

env_label={'Familiar','Novel','Familiar2'};
npil=[];
nothing=[];
missing=[];
% nothing=[15 20 25]; %E46
% nothing=[4 9 10 15 22]; %X14
% missing=[9]; %E46 days 4 5
% missing=[11]; %x14 days 2 3 4
% % Fcell{planes}=0;
Fall{3}=[];
forward{3}=[];
rotation{3}=[];
ybin{3}=[];
times{3}=[];
pos{3}=[];
% Features{planes}=0;
% Probs{planes}=0;

% v=[];
%%

for p=1:planes
    if manual==1;
        disp(['Plane ', num2str(p)]);
        [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    else
        names{p}=varargin{4}{p};
        paths{p}=varargin{5}{p};
    end
    %     files{p}=name{p}(1:end-4);
    load([paths{p},names{p}]);
%     lp=designfilt('lowpassfir','FilterOrder',128, ...
%         'PassbandFrequency',5,'StopbandFrequency',7,...
%         'SampleRate',length(forwardvel)/length(F)*Fs);
    F=Fc2;
    %     F=bsxfun(@rdivide,F,max(F,[],1));
    %     F=bsxfun(@minus,F,mean(F,1));
    %Downsample ybinned to match F
    if length(ybinned) ~= size(F,1)
        rewratio=length(ybinned)/size(F,1);
        ratio=floor(length(ybinned)/size(F,1));
        path=downsample(ybinned(1:(ratio*size(F,1))),ratio);
        
        fornan=forwardvel;
        fornan(abs(fornan)<.02)=NaN;
        rotnan=rotationvel;
        rotnan(abs(rotnan)<.02)=NaN;
        tidx=1:length(forwardvel);
        forwardvel=interp1(tidx(~isnan(fornan)),fornan(~isnan(fornan)),tidx,'pchip','extrap');
        rotationvel=interp1(tidx(~isnan(rotnan)),rotnan(~isnan(rotnan)),tidx,'pchip','extrap');
        for jj=1:length(F)
            if (jj*rewratio)<length(ybinned)
                rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
                forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
                ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            else
                rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):end));
                forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):end));
                ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):end));
            end
        end
        rotationvel((jj+1):end)=[];
        forwardvel((jj+1):end)=[];
        forwardvel=smooth(forwardvel);
        rotationvel=smooth(rotationvel);
        ybinned((jj+1):end)=[];
    end
    
    %Split into novel and familiar
    
    novel_start=round(timeSplit(1)/rewratio);
    novel_end=round(timeSplit(2)/rewratio);
    
    envinds{1}=1:novel_start;
    envinds{2}=(novel_start+1):novel_end;
    envinds{3}=(novel_end+1):length(F);
    for env=1:length(env_label)
        forward{env}=forwardvel(envinds{env});
        rotation{env}=rotationvel(envinds{env});
        ybin{env}=ybinned(envinds{env});
        paths{env}=path(envinds{env})-min(path);
        paths{env}=paths{env}/max(paths{env})*180/5;
        pos{env}=ceil(paths{env}+eps);
        times{env}=linspace(1,length(ybin{env})/Fs,length(ybin{env}));
        
        if ~isempty(Fall{env})
            Fall{env}=[Fall{env},F(envinds{env},:)];
        else
            Fall{env}=F(envinds{env},:);
        end
    end
end

if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
% Fall(del_cells)=[];

num_locs=max(pos{1});
num_cells=size(Fall{1},2);

%%
Falltem=[];
for env=1:length(env_label)
    dFcorr(num_cells,2)=0;
    for cell=1:num_cells
        [ctemp,ptemp]=corrcoef((Fall{env}(:,cell)),(forward{env}));
        dFcorr(cell,:)=[ctemp(2,1),ptemp(2,1)];
    end

%% Scatter dF and speed metrics
speedy=[0;diff(ybin{env})];
dFcorr(num_cells,2)=0;
for cell=1:num_cells
    [ctemp,ptemp]=corrcoef((Fall{env}(:,cell)),(forward{env}));
    dFcorr(cell,:)=[ctemp(2,1),ptemp(2,1)];
end

num_square=9;
num_figs=ceil(num_cells/num_square);
Falltem=bsxfun(@minus,Fall{env},min(Fall{env},[],1));
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
        scatter(log(abs(speedy)),log(Falltem(:,cell)),'b+')
        title(['F vs VR speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FVR scatter', num2str(f),env_label{env}, '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(forward{env})),log(Falltem(:,cell)),'b+')
        title(['F vs  Forward Ball Speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FBall scatter', num2str(f),env_label{env}, '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(forward{env})+abs(rotation{env})),log(Falltem(:,cell)),'b+')
        title(['F vs  Total Ball Speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FTotalBall scatter', num2str(f),env_label{env}, '.jpg']);
    
end
%test with 1 second bins
bin=round(Fs);
fbinwidth=(length(forward{env})-mod(length(forward{env}),bin))/bin;
binned = @(var) reshape(var(1:(fbinwidth*bin)),bin,fbinwidth);

favg=mean(binned(forward{env}))';
ravg=mean(binned(rotation{env}))';
vavg=mean(binned(speedy))';
fallavg=squeeze(mean(reshape(Fall{env}(1:(fbinwidth*bin),:),bin,fbinwidth,num_cells)));

num_square=9;
num_figs=ceil(num_cells/num_square);
% Falltem=bsxfun(@minus,Fall,min(Fall,[],1));
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
        scatter(log(abs(vavg)),(fallavg(:,cell)),'b+')
        title(['F vs VR speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FVR scatter binned ', num2str(f),env_label{env}, '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(favg)),(fallavg(:,cell)),'b+')
        title(['F vs  Forward Ball Speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FBall scatter binned ', num2str(f),env_label{env}, '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(sqrt((favg).^2+(ravg).^2)),(fallavg(:,cell)),'b+')
        title(['F vs  Total Ball Speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FTotalBall scatter binned ', num2str(f) ,env_label{env},'.jpg']);
    
end
    
    
    
%     figure
%     hold on
%     corrs=dFcorr(:,1);
%     ps=dFcorr(:,2);
%     
%     bar(corrs.*double(ps<.0005),'g')
%     bar(corrs.*double(ps>.0005),'r')
%     set(gca,'XTick',1:num_cells)
%     xlabel('Cell Number')
%     ylabel('Correlation')
%     title(['Forward Velocity and dF/F ',env_label{env}])
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'speed_dF',env_label{env}, '.jpg']);
%     
%     win=hanning(2^(floor(log2(length(Fall{env})))-1));
%     nover=2^(floor(log2(length(Fall{env})))-2);
%     nff=2^(floor(log2(length(Fall{env})))-1);
%     %% Analyze speed - use diff(ybinned) and sum(forv yawv) compare with activity
%     figure
%     speedy=[0;diff(ybin{env})];
%     speedv=forward{env};%+rotationvel;
%     plot(times{env},speedy/max(speedy),times{env},speedv/max(speedv)+1);
%     legend('VR Speed','Ball forward+yaw')
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR and Ball Speed',env_label{env},'.jpg']);
%     
%     figure
%     pwelch(speedy,win,nover,nff,Fs)
%     title('PSD Virtual reality speed')
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR PSD',env_label{env},'.jpg']);
%     
%     
%     
%     figure
%     pwelch(speedv,win,nover,nff,Fs)
%     title('PSD Ball speed')
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Ball PSD',env_label{env},'.jpg']);
%     
%     
%     
%     figure
%     pwelch(ybinned,win,nover,nff,Fs)
%     title('PSD Virtual reality position')
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Position PSD',env_label{env},'.jpg']);
%     
%     
%     %% xcov with speed
%     plot_xcorr1(Fall{env},speedv,Fs,saveDir,MouseID,[' dF and ball forward speed',env_label{env}]);
%     close all
%     [cross_corrs,~]=plot_xcorr1(Fall{env},sqrt(speedv.^2+rotation{env}.^2),Fs,saveDir,MouseID,[' dF and ball total speed',env_label{env}]);
%     plot_xcorr1(Fall{env},[0;diff(sqrt(speedv.^2+rotation{env}.^2))],Fs,saveDir,MouseID,[' dF and ball total accel',env_label{env}]);
%     close all
%     plot_xcorr1(Fall{env},speedy,Fs,saveDir,MouseID,[' dF and vr speed',env_label{env}]);
%     plot_xcorr1(Fall{env},[0;diff(speedy)],Fs,saveDir,MouseID,[' dF and vr accel',env_label{env}]);
%     close all
%     
%     %% xcorr cells
%     
%     plot_xcorr2(Fall{env},Fs,saveDir,[MouseID,' ',env_label{env}]);
%     close all
%     
%     %% power spectrum
%     [pwn,fpwn]=pwelch(Fall{env}(:,:),win,nover,nff,Fs);
%     
%     %% plotting
%     num_square=9;
%     num_figs=ceil(num_cells/num_square);
%     for f=1:num_figs
%         fig_start=(f-1)*num_square+1;
%         if f<num_figs
%             fig_end=fig_start+(num_square-1);
%         else fig_end = num_cells;
%         end
%         
%         %psd
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         suptitle(['PSD for cells: ',num2str(fig_start),' to ',num2str(fig_end), ' Mouse: ',MouseID])
%         
%         splot=1;
%         for cell=fig_start:fig_end
%             subplot(sqrt(num_square),sqrt(num_square),splot)
%             hold on
%             len=round(length(fpwn)/10);
%             plot(fpwn,(10*log10(pwn(:,cell))))
%             xlim([0 .05])
%             if ismember(cell,npil)
%                 title('Neurites')
%             elseif ismember(cell, nothing)
%                 title('Nothing Visible')
%             elseif ismember(cell,missing)
%                 title(['(Missing) Cell ', num2str(cell)])
%             else
%                 title(['Cell ', num2str(cell)])
%             end
%             splot=splot+1;
%             
%         end
%         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD_small', num2str(f),env_label{env}, '.jpg']);
%         savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD_small', num2str(f),env_label{env},'.fig']);
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         suptitle(['PSD for cells: ',num2str(fig_start),' to ',num2str(fig_end), ' Mouse: ',MouseID])
%         
%         splot=1;
%         for cell=fig_start:fig_end
%             subplot(sqrt(num_square),sqrt(num_square),splot)
%             hold on
%             len=round(length(fpwn)/10);
%             plot(fpwn,(10*log10(pwn(:,cell))))
%             if ismember(cell,npil)
%                 title('Neurites')
%             elseif ismember(cell, nothing)
%                 title('Nothing Visible')
%             elseif ismember(cell,missing)
%                 title(['(Missing) Cell ', num2str(cell)])
%             else
%                 title(['Cell ', num2str(cell)])
%             end
%             splot=splot+1;
%             
%         end
%         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f), env_label{env},'.jpg']);
%         savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f),env_label{env},'.fig']);
%         
%     end
%     
%     %% connection test
%     order=1:num_cells;
%     % order=[1 3 21 2 8 18 19 22 4 7 9 17 20 12 13];
%     % sorder=strsplit(num2str([order, npil, nothing]));
%     figure
%     % imagesc(corrcoef(Fall(:,[order, npil, nothing])))
%     title(['Correlation Matrix for mouse: ',MouseID])
%     imagesc(corrcoef(Fall{env}))
%     
%     colorbar
%     % set(gca,'Xtick',1:length(sorder),'XTickLabel',sorder)
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'CorrMat',env_label{env},'.jpg']);
%     
%     %% Filter and freq analysis
%     lsig0=zeros(size(Fall{env}));
%     hsig=lsig0;
%     gsig=hsig;
%     asig=gsig;
%     num_square=4;
%     num_figs=ceil(num_cells/num_square);
%     for f=1:num_figs
%         fig_start=(f-1)*num_square+1;
%         if f<num_figs
%             fig_end=fig_start+(num_square-1);
%         else fig_end = num_cells;
%         end
%         
%         %psd
%         %     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         splot=1;
%         for cell=fig_start:fig_end
%             
%             %         subplot(sqrt(num_square),sqrt(num_square),splot)
%             hold on
%             
%             lp=designfilt('lowpassfir','FilterOrder',128, ...
%                 'PassbandFrequency',.2,'StopbandFrequency',.4,...
%                 'SampleRate',Fs);
%             lsig0(:,cell)=filtfilt(lp,(Fall{env}(:,cell)));
%             
%             
%             hp=designfilt('highpassfir','FilterOrder',128, ...
%                 'PassbandFrequency',.2,'StopbandFrequency',.1,...
%                 'SampleRate',Fs);
%             if Fs==100
%                 gp=designfilt('highpassiir','FilterOrder',20, ...
%                     'PassbandFrequency',40,'PassbandRipple',.01,...
%                     'SampleRate',Fs);
%                 gsig(:,cell)=filtfilt(gp,Fall{env}(:,cell));
%                 ap=designfilt('bandpassiir','FilterOrder',20, ...
%                     'HalfPowerFrequency1',10,'HalfPowerFrequency2',20,...
%                     'PassbandRipple',.01,'SampleRate',Fs);
%                 asig(:,cell)=filtfilt(ap,Fall{env}(:,cell));
%             end
%             
%             hsig(:,cell)=filtfilt(hp,Fall{env}(:,cell));
%             %         plot(times,Fall(:,cell)+max(Fall(:,cell))*2,times,hsig+max(Fall(:,cell)),times,lsig0)
%             %         xlim([120*4 120*5])
%             %         legend('Cell Trace','High F','Low F')
%             %         plot(times,pos/max(pos)+(min(Fall(:,cell))-1))
%             
%             splot=splot+1;
%             
%             
%         end
%         %     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'freq bands', num2str(f), '.jpg']);
%         %     savefig([saveDir,'\',MouseID,'\', MouseID, 'freq bands', num2str(f),'.fig']);
%     end
%     
%     % [Coh,Fcoh]=mscohere(Fall(:,1),Fall(:,3),hann(512),128,512,Fs);
%     
%     
%     %% More freq analysis
%     params.fpass = [0 Fs/2];
%     params.Fs=Fs;
%     params.tapers=[5 9];
%     movingwin = [5 1];
%     for jj=1:num_cells
%         % jj=1;
%         %Spectrogram
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         suptitle(['Spectrogram for cell: ',num2str(jj), ' Mouse: ',MouseID])
%         
%         a1=subplot(4,1,1:2);
%         [specgram_data.s,specgram_data.t,specgram_data.f]=mtspecgramc(Fall{env}(:,jj),movingwin,params);
%         imagesc(specgram_data.t,specgram_data.f, log(specgram_data.s'))
%         
%         
%         a2=subplot(4,1,3:4);
%         hold on
%         plot(times{env},pos{env}/max(pos{env})+(max(Fall{env}(:,jj))))
%         
%         plot(times{env},speedv/max(speedv)+(min(Fall{env}(:,jj))-1))
%         plot(times{env},Fall{env}(:,jj));
%         xlim([0 max(times{env})])
%         linkaxes([a1,a2],'x');
%         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'spectrogram Tres', num2str(jj), env_label{env},'.jpg']);
%         savefig([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Tres', num2str(jj),env_label{env},'.fig']);
%     end
%     movingwin = [50 20];
%     for jj=1:num_cells
%         % jj=1;
%         %Spectrogram
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         suptitle(['Spectrogram for cell: ',num2str(jj), ' Mouse: ',MouseID])
%         
%         a1=subplot(4,1,1:2);
%         [specgram_data.s,specgram_data.t,specgram_data.f]=mtspecgramc(Fall{env}(:,jj),movingwin,params);
%         imagesc(specgram_data.t,specgram_data.f, log(specgram_data.s'))
%         
%         
%         a2=subplot(4,1,3:4);
%         hold on
%         plot(times{env},pos{env}/max(pos{env})+(max(Fall{env}(:,jj))))
%         
%         plot(times{env},speedv/max(speedv)+(min(Fall{env}(:,jj))-1))
%         plot(times{env},Fall{env}(:,jj));
%         xlim([0 max(times{env})])
%         linkaxes([a1,a2],'x');
%         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'spectrogram Fres', num2str(jj),env_label{env}, '.jpg']);
%         savefig([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Fres', num2str(jj),env_label{env},'.fig']);
%     end
%     
%     
%     %% Slow oscillations with speed, ybinned
%     swin=120;
%     ewin=240;
%     sframe=round(swin*Fs);
%     eframe=round(ewin*Fs);
%     constrain = @(sig) (sig-min(sig(sframe:eframe)))/max(sig(sframe:eframe)-min(sig(sframe:eframe)));
%     for cell=1:num_cells
%         gmag=abs(hilbert(gsig(:,cell)));
%         amag=abs(hilbert(asig(:,cell)));
%         
%         %     for n=2:2
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         %     offset=max(Fall(:,cell))/max(Fall(sframe:eframe,cell));
%         offset=1;
%         suptitle(['Oscillations for cell: ',num2str(cell), ' Mouse: ',MouseID])
%         a1=subplot(3,2,1);
%         hold on
%         %     plot(times{env},Fall(:,cell)/max(Fall(sframe:eframe,cell))+offset*2,...
%         %         times{env},hsig(:,cell)/max(hsig(sframe:eframe,cell))+offset,...
%         %         times{env},lsig0(:,cell)/max(lsig0(sframe:eframe,cell)));
%         %     if Fs==100
%         %         plot(times{env},gsig(:,cell)/max(gsig(sframe:eframe,cell))+3*offset);
%         %     end
%         %     plot(times{env},pos/max(pos)+min(Fall(sframe:eframe,1)-1));%,times{env},(Fall(:,1)))
%         plot(times{env},constrain(Fall{env}(:,cell))+offset*2,...
%             times{env},constrain(hsig(:,cell))+offset,...
%             times{env},constrain(lsig0(:,cell)));
%         if Fs==100
%             %         plot(times{env},constrain(gsig(:,cell))+3*offset);
%             plot(times{env},constrain(gmag)+3*offset);
%             plot(times{env},constrain(amag)+4*offset);
%         end
%         plot(times{env},pos{env}/max(pos{env})-1);%,times{env},(Fall(:,1)))
%         xlim([swin ewin])
%         title('Activity and Position')
%         a2=subplot(3,2,3);
%         hold on
%         plot(times{env},constrain(Fall{env}(:,cell))+offset*2,...
%             times{env},constrain(hsig(:,cell))+offset,...
%             times{env},constrain(lsig0(:,cell)));
%         if Fs==100
%             plot(times{env},constrain(gmag)+3*offset);
%             plot(times{env},constrain(amag)+4*offset);
%         end
%         plot(times{env},speedv/max(speedv)-1);%,times{env},(Fall(:,1)))
%         title('Activity and VR speed')
%         xlim([swin ewin])
%         
%         a3=subplot(3,2,5);
%         hold on
%         plot(times{env},constrain(Fall{env}(:,cell))+offset*2,...
%             times{env},constrain(hsig(:,cell))+offset,...
%             times{env},constrain(lsig0(:,cell)));
%         if Fs==100
%             plot(times{env},constrain(gmag)+3*offset);
%             plot(times{env},constrain(amag)+4*offset);
%         end
%         plot(times{env},speedy/max(speedy)-1);%,times{env},(Fall(:,1)))
%         title('Activity and Ball Speed')
%         xlim([swin ewin])
%         
%         subplot(3,2,2)
%         plot(fpwn,(10*log10(pwn(:,cell))))
%         
%         subplot(3,2,4)
%         plot(fpwn,(10*log10(pwn(:,cell))))
%         xlim([0 .05])
%         
%         
%         subplot(3,2,6)
%         plot(cross_corrs(round(size(cross_corrs,1)/2-Fs*2):round(size(cross_corrs,1)/2+Fs*2),cell))
%         title('Cross Corr between cell and ball speed')
%         linkaxes([a1,a2,a3],'x');
%         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'close look ', num2str(cell),',',env_label{env}, '.jpg']);
%         savefig([saveDir,'\',MouseID,'\', MouseID, 'close look', num2str(cell),',',env_label{env},'.fig']);
%         %     end
%     end
%     
%     env_label{env}
%     [feature_count,zM]=find_places(Fall{env},ybin{env},MouseID);
%     if isstruct(feature_count)
%     plot_places(feature_count.up,feature_count.down, saveDir, MouseID, env_label{env});
%     plot_places_MI(feature_count.up,zM.up(1,:),saveDir,MouseID,[env_label{env},' up'])
%     plot_places_MI(feature_count.down,zM.down(1,:),saveDir,MouseID,[env_label{env},' down'])
%     end
%     %%
%     close all
%     Fallt=[Fallt;Fall{env}];
end
% calc_PRTA(Fallt,rewards,round(Fs),MouseID,saveDir);
