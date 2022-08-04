%% Calculates PRTA(Reward synced activity) input=Fall,Rewards, Fs, MouseID, saveDir
% output=PRTA,meanPRTA,semPRTA

function [cPRTA,cmeanPRTA,csemPRTA]=calc_PRTA(Fall,rewards,Fs,MouseID,saveDir)
if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
ratio=floor(length(rewards)/size(Fall,1));
rewsmall2=zeros(length(Fall),1);
num_cells=size(Fall,2);
for jj=1:length(Fall)
    rewsmall2(jj)=max(rewards(((jj-1)*ratio+1):(ratio*jj)));
end
[~,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs)*2);
avgrewdist=median(diff(rewlocs));

timelook=[1,round(avgrewdist/Fs*2/3),16];
cPRTA{length(timelook)}=0;
cmeanPRTA{length(timelook)}=0;
csemPRTA{length(timelook)}=0;
for tl=1:length(timelook)
Fspan=round(Fs)*timelook(tl);
% rewlocs=find(rews);
cPRTA{tl} = zeros(Fspan*4+1,length(rewlocs),num_cells);
cmeanPRTA{tl}=zeros(Fspan*4+1,num_cells);
csemPRTA{tl}=cmeanPRTA{tl};
% Fall(:,end+1)=zscore(squeeze(mean(mean(video,1),2)));
% num_cells=num_cells+1;
for cell=1:num_cells
    for r=1:length(rewlocs)
        if rewlocs(r)<=Fspan*2
            cPRTA{tl}(:,r,cell)=[zeros(Fspan*2-rewlocs(r)+1,1);Fall(1:(rewlocs(r)+Fspan*2),cell)];
        elseif rewlocs(r)>=(length(Fall)-Fspan*2)
%             PRTA(:,r,cell)=[Fall((rewlocs(r)-Fspan*2):end,cell);zeros(rewlocs(r)+Fspan*2-length(rewsmall),1)];
                cPRTA{tl}(:,r,cell)=[Fall((rewlocs(r)-Fspan*2):end,cell);zeros(size(cPRTA{tl},1)-length(Fall((rewlocs(r)-Fspan*2):end,cell)),1)];
        else
            cPRTA{tl}(:,r,cell)=Fall((rewlocs(r)-Fspan*2):(rewlocs(r)+Fspan*2),cell);
        end
    end
    cmeanPRTA{tl}(:,cell)=mean(squeeze(cPRTA{tl}(:,:,cell)),2);
    csemPRTA{tl}(:,cell)=std(squeeze(cPRTA{tl}(:,:,cell)),0,2)/sqrt(size(squeeze(cPRTA{tl}(:,:,cell)),2));
end
end




for c=1:length(cPRTA)
    PRTA=cPRTA{c};
    meanPRTA=cmeanPRTA{c};
    semPRTA=csemPRTA{c};
    num_square=9;
    num_figs=ceil(num_cells/num_square);
    tl=(size(PRTA,1)-1)/(4*Fs);
    Fspan=round(Fs)*tl;
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
        imagesc(squeeze(PRTA(:,:,cell))')
        xlim([1 size(PRTA,1)])
        set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
        xlabel('Seconds')
        ylim([1 size(PRTA,2)])
        title(['PRTA heatmap cell: ', num2str(cell),' mouse ', MouseID]);
        splot=splot+1;
        
    end
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap', num2str(f), ' ', num2str(tl*4),' seconds','.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap', num2str(f), ' ', num2str(tl*4),' seconds','.fig']);

    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        errorbar(meanPRTA(:,cell),semPRTA(:,cell))
        line([2*Fspan,2*Fspan],[min(meanPRTA(:,cell)-3*semPRTA(:,cell)),max(meanPRTA(:,cell)+2*semPRTA(:,cell))]);
        
        xlim([0 size(meanPRTA,1)])
        set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
        xlabel('Seconds')
        title(['PRTA errorbar cell: ', num2str(cell),' mouse ', MouseID]);
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA errorbar', num2str(f), ' ', num2str(tl*4),' seconds', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA errorbar', num2str(f), ' ', num2str(tl*4),' seconds','.fig']);
    
%     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        hold on
        for nr=1:length(rewlocs)
            p=plot(PRTA(:,nr,cell)+(nr-1)*max(PRTA(:)));
            otherrews=intersect(find(rewlocs<(rewlocs(nr)+Fspan*2)),find(rewlocs>(rewlocs(nr)-Fspan*2)));
            s=scatter(rewlocs(otherrews)+2*Fspan-rewlocs(nr),Fall(rewlocs(otherrews)-1,cell)+(nr-1)*max(PRTA(:)),144,'+','Linewidth',3);
%             set(s,'CData',p.Color);
        end
%             errorbar(meanPRTA(:,cell),semPRTA(:,cell))
        line([2*Fspan,2*Fspan],[0,nr*max(PRTA(:))]);
        
        xlim([0 size(meanPRTA,1)])
        set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
        xlabel('Seconds')
        title(['PRTA traces cell: ', num2str(cell),' mouse ', MouseID]);
        splot=splot+1;
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA all traces cell ', num2str(cell), ' ', num2str(tl*4),' seconds', '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA all traces', num2str(f), ' ', num2str(tl*4),' seconds','.fig']);
    
    

    end
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    %plot all F with lines to show split
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
%         if Fs==15.5/4
%             t=0:1/Fs:(4*Fspan/Fs);%+1/Fs;
%         else
%             t=0:1/Fs:(4*Fspan/Fs);
%         end
        t=0:(size(meanPRTA,1)-1);
        plot(t,meanPRTA(:,cell),'b',t,(meanPRTA(:,cell)+2*semPRTA(:,cell)),'g:',t,(meanPRTA(:,cell)-2*semPRTA(:,cell)),'g:')
        set(gca,'XTick',[t(1:Fspan:end)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
        xlabel('Seconds')
        line([2*Fspan,2*Fspan],[min(meanPRTA(:,cell)-3*semPRTA(:,cell)),max(meanPRTA(:,cell)+2*semPRTA(:,cell))]);
        line([2*Fspan-avgrewdist,2*Fspan-avgrewdist],[min(meanPRTA(:,cell)-3*semPRTA(:,cell)),max(meanPRTA(:,cell)+2*semPRTA(:,cell))],'Color','r');
        line([2*Fspan+avgrewdist,2*Fspan+avgrewdist],[min(meanPRTA(:,cell)-3*semPRTA(:,cell)),max(meanPRTA(:,cell)+2*semPRTA(:,cell))],'Color','r');
        xlim([0 size(PRTA,1)])

        title(['PRTA confidence cell: ', num2str(cell),' mouse ',MouseID]);
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA confidence', num2str(f), ' ', num2str(tl*4),' seconds', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA confidence', num2str(f), ' ', num2str(tl*4),' seconds','.fig']);
    
    

end
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        [~,prtainds]=max(meanPRTA);
        [~,sortedinds]=sort(prtainds);
        imagesc(meanPRTA(:,sortedinds)')
        set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
        xlabel('Seconds')

        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap all', ' ', num2str(tl*4),' seconds', num2str(f), '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap all', ' ', num2str(tl*4),' seconds', num2str(f),'.fig']);
    
    
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         surf(meanPRTA(:,sortedinds)')
    
% num_square=ceil(sqrt(num_cells))^2;
% num_figs=ceil(num_cells/num_square);
% for cell=1:num_cells
%     rtemp=zeros(8*Fspan+1,sum(rewsmall));
%     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%     splot=1;
%     
%     for n = 1:num_cells
%         for r=1:sum(rewsmall)
%             rtemp(:,r)=xcov(PRTA(:,r,cell),PRTA(:,r,n),'coeff');
%         end
%         
%         subplot(sqrt(num_square),sqrt(num_square),splot)
%         imagesc(rtemp)
%         title(['PRTA xcov cell: ', num2str(cell),' and cell: ',num2str(n),' Mouse ', MouseID]);
%         splot=splot+1;
%         
%     end
%     
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA Xcorr cell', ' ', num2str(tl*4),' seconds', num2str(cell), '.jpg']);
%     savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA Xcorr cell', ' ', num2str(tl*4),' seconds', num2str(cell),'.fig']);
% end
end