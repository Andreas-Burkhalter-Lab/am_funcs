%Figures needed:
% SST
% % % 1 trace of sst behaviors aligned with behavior (movement?)
% 1 total plane F
% picture of selected cells (with numbers by them)
% % % df of selected cells in epochs, each day
% avg of all neurons, each day
% similar for more SSTs
% 
% AVG of ALL SSTs each day
% 
% PV
% entire plane dF
% individual cells D1-5
% identify good pv movie
% 

clear all;
%Add FOV df



%% Load and declare
% function interneuron_tests(MouseID, planes, paths, name)
saveDir='F:\MA Data\Interneurons\Figs\';
Fs=100;
planes=input('Number of planes: ');
days=input('Number of days: ');
% planes=4;
MouseID = input('Mouse ID: ');
% MouseID='E31D1';
names{days,planes}=0;
% names={'1119_000_000_hmm1_roibyhand_F_within.mat','1119_000_000_hmm2_roibyhand_F_within.mat',...
%     '1119_000_000_hmm3_roibyhand_F_within.mat','1119_000_000_hmm4_roibyhand_F_within.mat'};
paths{days,planes}=0;
% paths={'E:\2015\E31\151119\','E:\2015\E31\151119\',...
%     'E:\2015\E31\151119\','E:\2015\E31\151119\'};
% num_files=1;

npil=[];
% nothing=[15 20 25]; %E46
% nothing=[4 9 10 15 22]; %X14
% missing=[9]; %E46 days 4 5
% missing=[11]; %x14 days 2 3 4 
missing=[];
% % Fcell{planes}=0;
Fall{days}=[];
% Features{planes}=0;
% Probs{planes}=0;

% v=[];
for d=1:days
for p=1:planes
    disp(['Day: ', num2str(d), ' Plane ', num2str(p)]);
    [names{d,p},paths{d,p}]=uigetfile('*.mat','pick your files');
    load([paths{d,p},names{d,p}]);
    F=Fc2;
    %Downsample ybinned to match F
    if length(ybinned) ~= size(F,1)
        rewratiodays{d}=length(ybinned)/size(F,1);
        ratiodays{d}=floor(length(ybinned)/size(F,1));
        pathdays{d}=downsample(ybinned(1:(ratiodays{d}*size(F,1))),ratiodays{d});
        forwardveldays{d}=downsample(forwardvel(1:(ratiodays{d}*size(F,1))),ratiodays{d});
        rotationveldays{d}=downsample(rotationvel(1:(ratiodays{d}*size(F,1))),ratiodays{d});
        
    end
    rewdays{d}=rewards;
    pathdays{d}=(pathdays{d}-min(pathdays{d}));
    pathdays{d}=pathdays{d}/max(pathdays{d})*180/5;
    posdays{d}=ceil(pathdays{d}+eps);
    if ~isempty(Fall{d})
        if length(Fall{d})<length(F)
            Fall{d}=[Fall{d},F(1:length(Fall{d}),:)];
        else
            Fall{d}=[Fall{d}(1:length(F),:), F];
        end
    else
        Fall{d}=F;
    end
    novel_start{d}=round(timeSplit(1)/rewratiodays{d});
    novel_end{d}=round(timeSplit(2)/rewratiodays{d});
end
end
if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
%%
%E31
% use=[1 3 6 10 11];
% del_cells{1}=[5 6 10 11 15 16 23 24];
% del_cells{2}=del_cells{1};
% del_cells{3}=del_cells{2};
% del_cells{4}=del_cells{3};
% del_cells{5}=del_cells{3};
%X14
use=[1 2 9 11 16];
del_cells{1}=[9 15 20 25 ];
del_cells{2}=[4 9 10 11 15 21];
del_cells{3}=del_cells{2};
del_cells{4}=del_cells{3};
del_cells{5}=del_cells{3};
for d=1:days
   Fall{d}(:,del_cells{d})=[];
end
% Fall(:,del_cells)=[];
num_locs=max(posdays{1});
num_cells=size(Fall{1},2);
save([saveDir,MouseID,'\','vars']);
%%

% SST
% 1 trace of sst behaviors aligned with behavior (movement?)
% 1 total plane F
% picture of selected cells (with numbers by them)
% df of selected cells in epochs, each day
% avg of all neurons, each day
% similar for more SSTs
% 
% AVG of ALL SSTs each day
% 
% PV
% entire plane dF
% individual cells D1-5
% identify good pv movie
% 

% clear all;
%Add FOV df


%trace of 1 cell with behavior

% figure
% hold on plot
% plot(posdays{1}(:)/max(posdays{1})+1)
% plot(Fall{1}(:,3))


% 
% figure
% hold on
% plot(posdays{1}(950:1500)/max(posdays{1})+1)
% line([0 (Fs*10)], [min(Fall{1}(:,11)) min(Fall{1}(:,11))],'Color', 'k','LineWidth',3); %10S
% plot(Fall{1}(950:1500,3))
% line([550 550], [-.5 .5],'Color', 'r','LineWidth',3) %1 dF/F
% set(gca,'xticklabel', [],'xtick',[],'yticklabel',[],'ytick',[]);

% Fs=Fs;
seconds=Fs;

tenseconds=strsplit(num2str((Fs*4):(Fs*4):(Fs*4*9)));
figure
hold on
plot(posdays{1}(950:1500)/max(posdays{1})+1)
% line([0 (Fs*10)], [min(Fall{1}(:,11)) min(Fall{1}(:,11))],'Color', 'k','LineWidth',3); %10S
plot(Fall{1}(950:1500,3))
line([550 550], [-.5 .5],'Color', 'r','LineWidth',3) %1 dF/F
set(gca,'xticklabel', tenseconds,'xtick',[(Fs*4)*(1:9)],'yticklabel',[],'ytick',[]);
xlabel('Seconds')
title('Slow Oscillation in Cell Fluorescence')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'OscillationswithBeh.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'OscillationswithBeh.fig']);

% df of selected cells in epochs, each day

% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% for d=1:days
%     for u=1:length(use)
%         subplot(days,length(use),u+(d-1)*length(use))
%         cellF=Fall{d}(:,use(u));
%         mcell(1)=mean(cellF(1:novel_start{d}));
%         mcell(2)=mean(cellF(novel_start{d}+1:novel_end{d}));
%         mcell(3)=mean(cellF((novel_end{d}+1):end));
%         semcell(1)=1.96*std(cellF(1:novel_start{d}))/sqrt(length(1:novel_start{d}));
%         semcell(2)=1.96*std(cellF(novel_start{d}+1:novel_end{d}))/sqrt(length(novel_start{d}+1:novel_end{d}));
%         semcell(3)=1.96*std(cellF((novel_end{d}+1):end))/sqrt(length(cellF((novel_end{d}+1):end)));
%         hold on
%         plot(mcell,'k:','LineWidth',1)
%         errorbar(mcell,semcell,'r.','LineWidth',2)
%         ylim([-.05 .4])
%         set(gca,'xtick',[1 2 3],'xticklabel', {'Familiar', 'Novel','Familiar'})
%         ylabel('Mean dF/F');
%     end
% end
% saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'MeandFSelected6.svg']);
% savefig([saveDir,'\',MouseID,'\', MouseID, 'MeandFSelected6.fig']);
% 
xticks=linspace(0,length(Fall{1}),5);
xticklab=strsplit(num2str(xticks/(Fs)));
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% for d=1:days
%     for u=1:length(use)
%         subplot(days,length(use),u+(d-1)*length(use))
%         cellF=Fall{d}(:,use(u));
%         hold on
%         plot(cellF)
%         line([novel_start{d} novel_start{d}], [min(cellF) max(cellF)],'color','r')
%         line([novel_end{d} novel_end{d}], [min(cellF) max(cellF)],'color','r')
%         set(gca,'xticklabel', xticklab,'xtick',xticks,'yticklabel',[],'ytick',[]);
%         xlim([0 length(Fall{d})])
%         xlabel('Seconds')
%     end
% end
% 
% saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFtrace6.svg']);
% savefig([saveDir,'\',MouseID,'\', MouseID, 'dFtrace6.fig']);
use=use(1:5);

figure('units','normalized', 'Position', [.01 .05 .98 .87]);
for d=1:days
    for u=1:length(use)
        subplot(days,length(use),u+(d-1)*length(use))
        cellF=Fall{d}(:,use(u));
        mcell(1)=mean(cellF(1:novel_start{d}));
        mcell(2)=mean(cellF(novel_start{d}+1:novel_end{d}));
        mcell(3)=mean(cellF((novel_end{d}+1):end));
        semcell(1)=1.96*std(cellF(1:novel_start{d}))/sqrt(length(1:novel_start{d}));
        semcell(2)=1.96*std(cellF(novel_start{d}+1:novel_end{d}))/sqrt(length(novel_start{d}+1:novel_end{d}));
        semcell(3)=1.96*std(cellF((novel_end{d}+1):end))/sqrt(length(cellF((novel_end{d}+1):end)));
        hold on
        plot(mcell,'k:','LineWidth',1)
        errorbar(mcell,semcell,'r.','LineWidth',2)
        ylim([-.05 .4])
        set(gca,'xtick',[1 2 3],'xticklabel', {'Familiar', 'Novel','Familiar'})
        ylabel('Mean dF/F');
    end
end
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'MeandFSelected5.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'MeandFSelected5.fig']);

figure('units','normalized', 'Position', [.01 .05 .98 .87]);
for d=1:days
    for u=1:length(use)
        subplot(days,length(use),u+(d-1)*length(use))
        cellF=Fall{d}(:,use(u));
        hold on
        plot(cellF)
        line([novel_start{d} novel_start{d}], [min(cellF) max(cellF)],'color','r')
        line([novel_end{d} novel_end{d}], [min(cellF) max(cellF)],'color','r')
        set(gca,'xticklabel', xticklab,'xtick',xticks);
        xlim([0 length(Fall{d})])
        ylim([-1 2])
        xlabel('Seconds')
        ylabel('dF/F')
    end
end


saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFtrace5.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dFtrace5.fig']);
use=use(1:4);
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
for d=1:days
    for u=1:length(use)
        subplot(days,length(use),u+(d-1)*length(use))
        cellF=Fall{d}(:,use(u));
        mcell(1)=mean(cellF(1:novel_start{d}));
        mcell(2)=mean(cellF(novel_start{d}+1:novel_end{d}));
        mcell(3)=mean(cellF((novel_end{d}+1):end));
        semcell(1)=1.96*std(cellF(1:novel_start{d}))/sqrt(length(1:novel_start{d}));
        semcell(2)=1.96*std(cellF(novel_start{d}+1:novel_end{d}))/sqrt(length(novel_start{d}+1:novel_end{d}));
        semcell(3)=1.96*std(cellF((novel_end{d}+1):end))/sqrt(length(cellF((novel_end{d}+1):end)));
        hold on
        plot(mcell,'k:','LineWidth',1)
        errorbar(mcell,semcell,'r.','LineWidth',2)
        ylim([-.05 .4])
        set(gca,'xtick',[1 2 3],'xticklabel', {'Familiar', 'Novel','Familiar'})
        ylabel('Mean dF/F');
    end
end

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'MeandFSelected4.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'MeandFSelected4.fig']);
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
for d=1:days
    for u=1:length(use)
        subplot(days,length(use),u+(d-1)*length(use))
        cellF=Fall{d}(:,use(u));
        hold on
        plot(cellF)
        line([novel_start{d} novel_start{d}], [min(cellF) max(cellF)],'color','r')
        line([novel_end{d} novel_end{d}], [min(cellF) max(cellF)],'color','r')
        set(gca,'xticklabel', xticklab,'xtick',xticks);
        xlim([0 length(Fall{d})])
        ylim([-1 2])
        xlabel('Seconds')
        ylabel('dF/F')
    end
end


saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFtrace4.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dFtrace4.fig']);
%%Avg of all Neurons
figure('units','normalized', 'Position', [.01 .05 .24 .87]);
for d=1:days
    subplot(days,1,d)
    sumcells{d}(1)=mean(mean(Fall{d}(1:novel_start{d},:)));
    sumcells{d}(2)=mean(mean(Fall{d}((novel_start{d}+1):novel_end{d},:)));
    sumcells{d}(3)=mean(mean(Fall{d}((novel_end{d}+1):end,:)));
    semcells{d}(1)=1.96*std(mean(Fall{d}(1:novel_start{d},:),2))/sqrt(length(1:novel_start{d}));
    semcells{d}(2)=1.96*std(mean(Fall{d}((novel_start{d}+1):novel_end{d},:)))/sqrt(length(novel_start{d}+1:novel_end{d}));
    semcells{d}(3)=1.96*std(mean(Fall{d}((novel_end{d}+1):end,:)))/sqrt(length(Fall{2}((novel_end{d}+1):end,1)));
  

%     sumcells{d}=sumcells{d}+sumcells1{d};
%     semcell{d}=semcell{d}+semcell1{d};
    hold on    
    errorbar(sumcells{d}/2,semcells{d}/2,'.r','markersize',40)
    plot(sumcells{d}/2,'k:')
    set(gca,'xtick',[1 2 3],'xticklabel', {'Familiar', 'Novel','Familiar'})
    ylim([-.05 .07])
end

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFallCell.svg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'dFallCell.fig']);
    save([saveDir,MouseID,'means.mat'],'sumcells','semcells');
    
    %load full fov for df
    [names,paths]=uigetfile('*.mat','pick your files');
    load([paths,names]); 
    xticks=linspace(0,length(F),5);
    xticklab=strsplit(num2str(xticks/(Fs)));
    fov=mean(mean(video,1),2);
    figure
    hold on
    plot(squeeze(fov))
    line([novel_start{1} novel_start{1}], [min(fov) max(fov)],'color','r')
    line([novel_end{1} novel_end{1}], [min(fov) max(fov)],'color','r')
    set(gca,'xticklabel', xticklab,'xtick',xticks)
    ylim([min(fov)*.99 max(fov)*1.01])
    xlabel('Seconds')
    ylabel('dF/F')
    title('Mean Fluorescence of Plane')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'fullFOV.svg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'fullFOV.fig']);
    figure
    mfov(1)=mean(fov(1:novel_start{1}));
    mfov(2)=mean(fov((novel_start{1}+1):novel_end{1}));
    mfov(3)=mean(fov((novel_end{1}+1):end));
    hold on
    title('Mean Fluorescence of Plane')
    errorbar(mfov,[0,0,0],'.r','LineWidth',10)
    plot(mfov,'k:')        
    ylabel('dF/F')
    set(gca,'xtick',[1 2 3],'xticklabel', {'Familiar', 'Novel','Familiar'})
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'fullFOVmean.svg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'fullFOVmean.fig']);
    figure
    imshow(mat2gray(max(video,[],3)));
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'fullFOV.jpg']);

%%
% 
% % dFcorr(num_cells,2)=0;
% % for cell=1:num_cells
% %     [ctemp,ptemp]=corrcoef((Fall(:,cell)),(v(:,cell)));
% %     dFcorr(cell,:)=[ctemp(2,1),ptemp(2,1)];
% % end
% % cell=13;
% % scatter(v(:,cell),Fall(:,cell))
% %
% %
% % figure
% % hold on
% % corrs=dFcorr(:,1);
% % ps=dFcorr(:,2);
% %
% % bar(corrs.*double(ps<.0005),'g')
% % bar(corrs.*double(ps>.0005),'r')
% % set(gca,'XTick',1:num_cells)
% % xlabel('Cell Number')
% % ylabel('Correlation')
% % title('Forward Velocity and dF/F')
% % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'speed_dF', '.svg']);
% 
% 
% %Analyze speed - use diff(ybinned) and sum(forv yawv) compare with activity
% % 
% % speedy{1}=(diff(ybinned(1:timeSplit(1))));
% % speedy{2}=(diff(ybinned((timeSplit(1)+1):timeSplit(2))));
% % speedy{3}=(diff(ybinned((timeSplit(2)+1):end)));
% % figure
% % plot([speedy{1}', speedy{2}', speedy{3}']);
% % title('Virtual reality speed')
% % avgspeedy(3)=0;
% % avgspeedv(3)=0;
% % semspeedy(3)=0;
% % semspeedv(3)=0;
% % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR Speed','.svg']);
% % pwelch([speedy{1}', speedy{2}',speedy{3}'],900000,450000,900000,Fs*rewratio)
% % 
% % speedv{1}=(forwardvel(1:novel_start))+rotationvel(1:novel_start);
% % speedv{2}=forwardvel((novel_start+1):novel_end)+rotationvel((novel_start+1):novel_end);
% % speedv{3}=forwardvel((novel_end+1):end)+rotationvel((novel_end+1):end);
% % figure
% % plot([speedv{1}', speedv{2}', speedv{3}']);
% % title('Ball speed')
% % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Ball Speed','.svg']);
% % 
% % for i=1:3
% %     avgspeedy(i)=sum(abs(speedy{i}))/length(speedy{i});
% %     semspeedy(i)=1.96*std(speedy{i})/sqrt(length(speedy{i}));
% %     avgspeedv(i)=sum(speedv{i})/length(speedv{i});
% %     semspeedv(i)=1.96*std(speedv{i})/sqrt(length(speedv{i}));
% % end
% % figure; hold on;
% % % bar(avgspeedy,'w')
% % errorbar(avgspeedy,semspeedy,'-x')
% % title('Average VR Speeds')
% % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Avg VR speeds','.svg']);
% % figure
% % errorbar(avgspeedv,semspeedv,'-x')
% % title('Average Ball Speeds')
% % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Avg Ball speeds','.svg']);
% 
% for d=1:days
% %Look at rew/min
% rewenv(1)=sum(double(rewdays{d}(1:timeSplit(1))));
% rewenv(2)=sum(double(rewdays{d}((timeSplit(1)+1):timeSplit(2))));
% rewenv(3)=sum(double(rewdays{d}((timeSplit(2)+1):end)));
% 
% rewenvminmax = max(rewenv(1)/(novel_start{d}/(Fs)/60));
% rewenvmin(1)=rewenv(1)/(novel_start{d}/(Fs)/60)/rewenvminmax;
% rewenvmin(2)=rewenv(2)/((novel_end{d}-novel_start{d})/(Fs)/60)/rewenvminmax;
% rewenvmin(3)=rewenv(3)/((length(F)-novel_end{d})/(Fs)/60)/rewenvminmax;
% 
% semrew(1)=1.96*std(rewenv(1))/sqrt(length(rewenv(1)));
% semrew(2)=1.96*std(rewenv(2))/sqrt(length(rewenv(2)));
% semrew(3)=1.96*std(rewenv(3))/sqrt(length(rewenv(3)));
% subplot(1,days+1,d)
% errorbar(rewenvmin,semrew)
% ylim([0 1.2])
% title('Rew/min')
% % set(gca,'XTickLabels',{'Familiar', 'Novel','Familiar'})
% novrew(d)=rewenvmin(2);
% novsem(d)=semrew(2);
% end
% subplot(1,days+1,days+1)
% errorbar(novrew,novsem)
% title('Rewards in Novel by Day')
% ylim([0 1.2])
% saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Rewpermins','.svg']);
% 
% % num_square=9;
% % mcell(3)=0;
% % semcell=mcell;
% % num_figs = ceil(num_cells/num_square);
% % %power spectrum
% % [pw,fpw]=pwelch(Fall,3000,1500,3000,Fs);
% % figure;
% % [poswelch,fposwelch]=pwelch(pos,3000,1500,3000,Fs);
% % plot(fposwelch,10*log10(poswelch))
% % 
% % figure;
% % [ywelch,fywelch]=pwelch(ybinned,900000,450000,900000,Fs*rewratio);
% % plot(fywelch(1:round(length(fywelch)/3000)),10*log10(ywelch(1:round(length(fywelch)/3000))))
% % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PositionSpectrum','.svg']);
% % 
% % %PSD for each environment
% % [pwn1,fpwn1]=pwelch(Fall(1:novel_start,:),1500,750,1500,Fs);
% % [pwn2,fpwn2]=pwelch(Fall(novel_start+1:novel_end,:),2500,1500,2500,Fs);
% % [pwn3,fpwn3]=pwelch(Fall(novel_end+1:end,:),1500,750,1500,Fs);
% 
% %%
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% for d=1:days
%     for cell=1:size(Fall{d},2)
%         
%         %         %plot means with sem
%         subplot(days,size(Fall{1},2),cell+(d-1)*size(Fall{2},2))    
%         hold on
%         cellF=Fall{d}(:,cell);
%         
% 
%         mcell(1)=mean(cellF(1:novel_start{d}));
%         mcell(2)=mean(cellF(novel_start{d}+1:novel_end{d}));
%         mcell(3)=mean(cellF((novel_end{d}+1):end));
%         
%         semcell(1)=1.96*std(cellF(1:novel_start{d}))/sqrt(length(1:novel_start{d}));
%         semcell(2)=1.96*std(cellF(novel_start{d}+1:novel_end{d}))/sqrt(length(novel_start{d}+1:novel_end{d}));
%         semcell(3)=1.96*std(cellF((novel_end{d}+1):end))/sqrt(length(cellF((novel_end{d}+1):end)));
%         plot(mcell,'k:')
%         errorbar(mcell,semcell,'r.')
%         ylim([-.1 .5])
% %         if ismember(cell,npil)
% %             title('Neurites')
% %         elseif ismember(cell, nothing)
% %             title('Nothing Visible')
% %         elseif ismember(cell,missing)
% %             title(['(Missing) Cell ', num2str(cell)])
% %         else
% %             title(['Cell ', num2str(cell)])
% %         end
% %         
%     end
% end
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenvMean.svg']);
%     savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenvBarMean.fig']);
%     
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% for d=1:days
%     for cell=1:size(Fall{d},2)
%         subplot(days,size(Fall{1},2),cell+(d-1)*size(Fall{2},2))    
%         hold on
%         cellF=Fall{d}(:,cell);
%         plot(cellF)
%         line([novel_start{d} novel_start{d}], [min(cellF) max(cellF)],'color','r')
%         line([novel_end{d} novel_end{d}], [min(cellF) max(cellF)],'color','r')
% %         if ismember(cell,npil)
% %             title('Neurites')
% %         elseif ismember(cell, nothing)
% %             title('Nothing Visible')
% %         elseif ismember(cell,missing)
% %             title(['(Missing) Cell ', num2str(cell)])
% %         else
% %             title(['Cell ', num2str(cell)])
% %         end
%         
%     end
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenv.svg']);
%     savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenv.fig']);
%     
%     %     for cell=fig_start:fig_end
%     %         subplot(sqrt(num_square),sqrt(num_square),splot)
%     %         hold on
%     %         cellF=Fall(:,cell);
%     %         nfft=2^nextpow2(length(cellF));
%     %         cellfft=fft(cellF,nfft)/nfft;
%     %         freq=(Fs*linspace(0,1,nfft/2+1));
%     %         plot(freq,2*abs(cellfft(1:nfft/2+1)))
%     %         if ismember(cell,npil)
%     %             title('Neurites')
%     %         elseif ismember(cell, nothing)
%     %             title('Nothing Visible')
%     %         else
%     %             title(['Cell ', num2str(cell)])
%     %         end
%     %         splot=splot+1;
%     %     end
%     %          saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f), '.svg']);
%     %          savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f)]);
% end
% 
% %% connection test
% order=1:size(Fall,2);
% % order=[1 3 21 2 8 18 19 22 4 7 9 17 20 12 13];
% % sorder=strsplit(num2str([order, npil, nothing]));
% 
% % imagesc(corrcoef(Fall(:,[order, npil, nothing])))
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% for d=1:days
%         subplot(1,days,d)    
%         imagesc(corrcoef(Fall{d}))
% end
% colormap('jet')
% % colorbar
% % set(gca,'Xtick',1:length(sorder),'XTickLabel',sorder)
% saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'CorrMat', '.svg']);
% 
% % 
% % lp=designfilt('lowpassiir','FilterOrder',8, ...
% %     'PassbandFrequency',.05,...
% %     'SampleRate',Fs);
% % lsig=filter(lp,temp);
% % lsig0=filtfilt(lp,temp);
% % 
% % figure
% % plot(path/max(path),lsig0)
% % figure
% % plot([path/max(path),lsig])
% % 
% % hp=designfilt('highpassiir','FilterOrder',8, ...
% %     'PassbandFrequency',.05,...
% %     'SampleRate',Fs);
% % hsig=filtfilt(hp,temp);
% % plot([temp,hsig,lsig])
% % 
% % henv=abs(hilbert(hsig));
% % figure; plot([path/max(path), henv])
% % 
% % 
% % [Coh,Fcoh]=mscohere(Fall(:,1),Fall(:,3),hann(512),128,512,Fs);
% % 
% 

