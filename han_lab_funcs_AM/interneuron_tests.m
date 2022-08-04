%Need to test activity vs velocity and activity vs environment

%will have 3 F-files for each plane

%cp code here to go through each

%code here for place cells

%here for df/f in f vs n
%make sure cells have the same index

%play mean F in each loc, tile 16 bars


%% Load and declare
clear all;
% function interneuron_tests(MouseID, planes, paths, name)
saveDir='G:\MA\Interneurons\';
Fs=15.5/4;
planes=input('Number of planes: ');
% planes=4;
MouseID = input('Mouse ID and Day: ');
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
    for p=1:planes
    disp(['Plane ', num2str(p)]);
    [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    %     files{p}=name{p}(1:end-4);
    load([paths{p},names{p}]);
    F=Fc2;
    %Downsample ybinned to match F
    if length(ybinned) ~= size(F,1)
        rewratio=length(ybinned)/size(F,1);
        ratio=floor(length(ybinned)/size(F,1));
        path=downsample(ybinned(1:(ratio*size(F,1))),ratio);
        forwardvel=downsample(forwardvel(1:(ratio*size(F,1))),ratio);
        rotationvel=downsample(rotationvel(1:(ratio*size(F,1))),ratio);
        
    end
    path=(path-min(path));
    path=path/max(path)*180/5;
    pos=ceil(path+eps);
    if ~isempty(Fall)
        if length(Fall)<length(F)
            Fall=[Fall,F(1:length(Fall),:)];
            v=[v, repmat(forwardvel(1:length(v)),1,size(F,2))];
        else
            Fall=[Fall(1:length(F),:), F];
            v=[v(1:length(forwardvel),:), repmat(forwardvel,1,size(F,2))];
        end
    else
        Fall=F;
        v=repmat(forwardvel,1,size(F,2));
    end
    novel_start=round(timeSplit(1)/rewratio);
    novel_end=round(timeSplit(2)/rewratio);
    end
    
    if ~exist([saveDir,MouseID,'\'],'dir')
        mkdir([saveDir,MouseID,'\']);
    end
    % Fall(del_cells)=[];
    num_locs=max(pos);
    num_cells=size(Fall,2);
    %%
    
    % dFcorr(num_cells,2)=0;
    % for cell=1:num_cells
    %     [ctemp,ptemp]=corrcoef((Fall(:,cell)),(v(:,cell)));
    %     dFcorr(cell,:)=[ctemp(2,1),ptemp(2,1)];
    % end
    % cell=13;
    % scatter(v(:,cell),Fall(:,cell))
    %
    %
    % figure
    % hold on
    % corrs=dFcorr(:,1);
    % ps=dFcorr(:,2);
    %
    % bar(corrs.*double(ps<.0005),'g')
    % bar(corrs.*double(ps>.0005),'r')
    % set(gca,'XTick',1:num_cells)
    % xlabel('Cell Number')
    % ylabel('Correlation')
    % title('Forward Velocity and dF/F')
    % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'speed_dF', '.jpg']);
    
    
    %% Analyze speed - use diff(ybinned) and sum(forv yawv) compare with activity
    
    speedy{1}=(diff(ybinned(1:timeSplit(1))));
    speedy{2}=(diff(ybinned((timeSplit(1)+1):timeSplit(2))));
    speedy{3}=(diff(ybinned((timeSplit(2)+1):end)));
    figure
    plot([speedy{1}', speedy{2}', speedy{3}']);
    title('Virtual reality speed')
    avgspeedy(3)=0;
    avgspeedv(3)=0;
    semspeedy(3)=0;
    semspeedv(3)=0;
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR Speed','.jpg']);
    pwelch([speedy{1}', speedy{2}',speedy{3}'],900000,450000,900000,Fs*rewratio)
    
    speedv{1}=(forwardvel(1:novel_start))+rotationvel(1:novel_start);
    speedv{2}=forwardvel((novel_start+1):novel_end)+rotationvel((novel_start+1):novel_end);
    speedv{3}=forwardvel((novel_end+1):end)+rotationvel((novel_end+1):end);
    figure
    plot([speedv{1}', speedv{2}', speedv{3}']);
    title('Ball speed')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Ball Speed','.jpg']);
    
    for i=1:3
        avgspeedy(i)=sum(abs(speedy{i}))/length(speedy{i});
        semspeedy(i)=1.96*std(speedy{i})/sqrt(length(speedy{i}));
        avgspeedv(i)=sum(speedv{i})/length(speedv{i});
        semspeedv(i)=1.96*std(speedv{i})/sqrt(length(speedv{i}));
    end
    figure; hold on;
    % bar(avgspeedy,'w')
    errorbar(avgspeedy,semspeedy,'-x')
    ylim([0 max(avgspeedy)*1.2])
    title('Average VR Speeds')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Avg VR speeds','.jpg']);
    figure
    errorbar(avgspeedv,semspeedv,'-x')
    ylim([0 max(avgspeedv)*1.2])
    title('Average Ball Speeds')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Avg Ball speeds','.jpg']);
    
    %% Look at rew/min
    rewenv(1)=sum(double(rewards(1:timeSplit(1))));
    rewenv(2)=sum(double(rewards((timeSplit(1)+1):timeSplit(2))));
    rewenv(3)=sum(double(rewards((timeSplit(2)+1):end)));
    
    rewenvmin(1)=rewenv(1)/(novel_start/(Fs)/60);
    rewenvmin(2)=rewenv(2)/((novel_end-novel_start)/(Fs)/60);
    rewenvmin(3)=rewenv(3)/((length(F)-novel_end)/(Fs)/60);
    
    semrew(1)=1.96*std(rewenv(1))/sqrt(length(rewenv(1)));
    semrew(2)=1.96*std(rewenv(2))/sqrt(length(rewenv(2)));
    semrew(3)=1.96*std(rewenv(3))/sqrt(length(rewenv(3)));
    errorbar(rewenvmin,semrew)
    ylim([0 max(rewenvmin)*1.2])
    title('Rew/min')
    % set(gca,'XTicks',[1 2 3],'XTickLabels',{'Familiar', 'Novel','Familiar'})
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Rewpermins','.jpg']);
    
    num_square=9;
    mcell(3)=0;
    semcell=mcell;
    num_figs = ceil(num_cells/num_square);
    %% power spectrum
    [pw,fpw]=pwelch(Fall,3000,1500,3000,Fs);
    figure;
    [poswelch,fposwelch]=pwelch(pos,3000,1500,3000,Fs);
    plot(fposwelch,(poswelch))
    
    figure;
    [ywelch,fywelch]=pwelch(ybinned,900000,450000,900000,Fs*rewratio);
    plot(fywelch(1:round(length(fywelch)/3000)),10*log10(ywelch(1:round(length(fywelch)/3000))))
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PositionSpectrum','.jpg']);
    
    %PSD for each environment
    [pwn1,fpwn1]=pwelch(Fall(1:novel_start,:),1500,750,1500,Fs);
    [pwn2,fpwn2]=pwelch(Fall(novel_start+1:novel_end,:),2500,1500,2500,Fs);
    [pwn3,fpwn3]=pwelch(Fall(novel_end+1:end,:),1500,750,1500,Fs);
    
    %% plotting
    
    for f=1:num_figs
        fig_start=(f-1)*num_square+1;
        if f<num_figs
            fig_end=fig_start+(num_square-1);
        else fig_end = num_cells;
        end
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        splot=1;
        %psd by parts
        for cell=fig_start:fig_end
            subplot(sqrt(num_square),sqrt(num_square),splot)
            hold on
            len1=round(length(fpwn1)/10);
            len2=round(length(fpwn2)/10);
            len3=round(length(fpwn3)/10);
            
            plot(fpwn1(1:len1),10*log10(pwn1(1:len1,cell)),'r')
            plot(fpwn2(1:len2),10*log10(pwn2(1:len2,cell)),'g')
            plot(fpwn3(1:len3),10*log10(pwn3(1:len3,cell)),'b')
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
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSDparts', num2str(f), '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'PSDparts', num2str(f),'.fig']);
        
        %psd
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        splot=1;
        for cell=fig_start:fig_end
            subplot(sqrt(num_square),sqrt(num_square),splot)
            hold on
            len=round(length(fpw)/10);
            plot(fpw(1:len),10*log10(pw(1:len,cell)))
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
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f), '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f),'.fig']);
        
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        splot=1;
        for cell=fig_start:fig_end
            %         %plot means with sem
            subplot(sqrt(num_square),sqrt(num_square),splot)
            hold on
            cellF=Fall(:,cell);
            
            
            mcell(1)=mean(cellF(1:novel_start));
            mcell(2)=mean(cellF(novel_start+1:novel_end));
            mcell(3)=mean(cellF((novel_end+1):end));
            
            semcell(1)=1.96*std(cellF(1:novel_start))/sqrt(length(1:novel_start));
            semcell(2)=1.96*std(cellF(novel_start+1:novel_end))/sqrt(length(novel_start+1:novel_end));
            semcell(3)=1.96*std(cellF((novel_end+1):end))/sqrt(length(cellF((novel_end+1):end)));
            plot(mcell,'k:')
            errorbar(mcell,semcell,'r.')
            ylim([-.1 .5])
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
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f), '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f),'.fig']);
        
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        splot=1;
        %plot all F with lines to show split
        for cell=fig_start:fig_end
            subplot(sqrt(num_square),sqrt(num_square),splot)
            hold on
            cellF=Fall(:,cell);
            plot(cellF)
            line([novel_start novel_start], [min(cellF) max(cellF)],'color','r')
            line([novel_end novel_end], [min(cellF) max(cellF)],'color','r')
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
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenv', num2str(f), '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenv', num2str(f),'.fig']);
        
        %     for cell=fig_start:fig_end
        %         subplot(sqrt(num_square),sqrt(num_square),splot)
        %         hold on
        %         cellF=Fall(:,cell);
        %         nfft=2^nextpow2(length(cellF));
        %         cellfft=fft(cellF,nfft)/nfft;
        %         freq=(Fs*linspace(0,1,nfft/2+1));
        %         plot(freq,2*abs(cellfft(1:nfft/2+1)))
        %         if ismember(cell,npil)
        %             title('Neurites')
        %         elseif ismember(cell, nothing)
        %             title('Nothing Visible')
        %         else
        %             title(['Cell ', num2str(cell)])
        %         end
        %         splot=splot+1;
        %     end
        %          saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f), '.jpg']);
        %          savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f)]);
    end
    
    %% connection test
    order=1:size(Fall,2);
    % order=[1 3 21 2 8 18 19 22 4 7 9 17 20 12 13];
    % sorder=strsplit(num2str([order, npil, nothing]));
    figure
    % imagesc(corrcoef(Fall(:,[order, npil, nothing])))
    imagesc(corrcoef(Fall))
    
    colormap('jet')
    colorbar
    % set(gca,'Xtick',1:length(sorder),'XTickLabel',sorder)
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'CorrMat', '.jpg']);
    
%     %% Filter and freq analysis
%     lp=designfilt('lowpassiir','FilterOrder',8, ...
%         'PassbandFrequency',.05,...
%         'SampleRate',Fs);
%     lsig=filter(lp,temp);
%     lsig0=filtfilt(lp,temp);
%     
%     figure
%     plot(path/max(path),lsig0)
%     figure
%     plot([path/max(path),lsig])
%     
%     hp=designfilt('highpassiir','FilterOrder',8, ...
%         'PassbandFrequency',.05,...
%         'SampleRate',Fs);
%     hsig=filtfilt(hp,temp);
%     plot([temp,hsig,lsig])
%     
%     henv=abs(hilbert(hsig));
%     figure; plot([path/max(path), henv])
%     
%     
%     [Coh,Fcoh]=mscohere(Fall(:,1),Fall(:,3),hann(512),128,512,Fs);
%     
%     
%     
