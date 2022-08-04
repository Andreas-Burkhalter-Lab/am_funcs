%% Calculates PRTposition(Reward synced position and speed and stuff) input=Fall,Rewards, Fs, MouseID, saveDir
% output=PRTA,meanPRTA,semPRTA

function [cPRTA,cmeanPRTA,csemPRTA,maxinfo]=calc_PRTloc(loc,rewards,Fs,MouseID,saveDir,varargin)
if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
if nargin>5
    env=varargin{1};
else env='';
end
if size(loc,2)>size(loc,1)
    loc=loc';
end
if nargin>6
    part=varargin{2};
end
rewratio=(length(rewards)/length(loc));
rewsmall2=zeros(length(loc),1);
% num_cells=size(loc,2);
for jj=1:length(loc)
    rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
end
[~,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs)*2);
avgrewdist=median(diff(rewlocs));

% timelook=[1,round(avgrewdist/Fs*2/3),16];
timelook=[2 4 16];
cPRTA{length(timelook)}=0;
cmeanPRTA{length(timelook)}=0;
csemPRTA{length(timelook)}=0;
for tl=1:length(timelook)
    Fspan=round(Fs)*timelook(tl);
    % rewlocs=find(rews);
    cPRTA{tl} = zeros(Fspan*4+1,length(rewlocs));
    cmeanPRTA{tl}=zeros(Fspan*4+1,1);
    csemPRTA{tl}=cmeanPRTA{tl};
    % Fall(:,end+1)=zscore(squeeze(mean(mean(video,1),2)));
    % num_cells=num_cells+1;
    for r=1:length(rewlocs)
        if rewlocs(r)<=Fspan*2
            cPRTA{tl}(:,r)=[zeros(Fspan*2-rewlocs(r)+1,1);loc(1:(rewlocs(r)+Fspan*2))];
        elseif rewlocs(r)>=(length(loc)-Fspan*2)
            %             PRTA(:,r,cell)=[Fall((rewlocs(r)-Fspan*2):end,cell);zeros(rewlocs(r)+Fspan*2-length(rewsmall),1)];
            cPRTA{tl}(:,r)=[loc((rewlocs(r)-Fspan*2):end);zeros(size(cPRTA{tl},1)-length(loc((rewlocs(r)-Fspan*2):end)),1)];
        else
            cPRTA{tl}(:,r)=loc((rewlocs(r)-Fspan*2):(rewlocs(r)+Fspan*2));
        end
    end
    cmeanPRTA{tl}(:)=nanmean(squeeze(cPRTA{tl}(:,:)),2);
    csemPRTA{tl}(:)=nanstd(squeeze(cPRTA{tl}(:,:)),0,2)/sqrt(size(squeeze(cPRTA{tl}(:,:)),2));
end


constrain = @(sig) (sig-min(sig))/(max(sig)-min(sig));


for c=1:length(cPRTA)
    PRTA=cPRTA{c};
    meanPRTA=cmeanPRTA{c};
    semPRTA=csemPRTA{c};
    %     num_square=9;
    %     num_figs=ceil(num_cells/num_square);
    tl=(size(PRTA,1)-1)/(4*Fs);
    Fspan=round(Fs)*round(tl);
    % for f=1:num_figs
    %     fig_start=(f-1)*num_square+1;
    %     if f<num_figs
    %         fig_end=fig_start+(num_square-1);
    %     else fig_end = num_cells;
    %     end
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    %     splot=1;
    %     for cell=fig_start:fig_end
    %         subplot(sqrt(num_square),sqrt(num_square),splot)
    hold on
    imagesc(squeeze(PRTA(:,:))')
    xlim([1 size(PRTA,1)])
    set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
    xlabel('Seconds')
    ylim([1 max(size(PRTA,2),2)])
    title(['PRTloc heatmap mouse ', MouseID]);
    
    %     end
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTloc heatmap', ' ', num2str(tl*4),' seconds',env,'.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTloc heatmap', ' ', num2str(tl*4),' seconds',env,'.fig']);
    
    %     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    %     splot=1;
    %     for cell=fig_start:fig_end
    %         subplot(sqrt(num_square),sqrt(num_square),splot)
    %         hold on
    %         errorbar(meanPRTA(:,cell),semPRTA(:,cell))
    %         line([2*Fspan,2*Fspan],[min(meanPRTA(:,cell)-3*semPRTA(:,cell)),max(meanPRTA(:,cell)+2*semPRTA(:,cell))]);
    %
    %         xlim([0 size(meanPRTA,1)])
    %         set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
    %         xlabel('Seconds')
    %         title(['PRTA errorbar cell: ', num2str(cell),' mouse ', MouseID]);
    %         splot=splot+1;
    %
    %     end
    %     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA errorbar', num2str(f), ' ', num2str(tl*4),' seconds', '.jpg']);
    %     savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA errorbar', num2str(f), ' ', num2str(tl*4),' seconds','.fig']);
    %
    
    %     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    % %individual traces
    %
    %         hold on
    %         for nr=1:length(rewlocs)
    %             trace=constrain(PRTA(:,nr));
    %             plot(trace+.5*(nr-1));
    %             otherrews=intersect(find(rewlocs<(rewlocs(nr)+Fspan*2)),find(rewlocs>(rewlocs(nr)-Fspan*2)));
    % %             rewlocs(otherrews)
    % %             rewlocs(nr)
    % %             Fspan*2
    %             scatter(rewlocs(otherrews)+2*Fspan-rewlocs(nr),trace(rewlocs(otherrews)-rewlocs(nr)+Fspan*2)+.5*(nr-1),144,'+','Linewidth',3);
    %
    % %             set(s,'CData',p.Color);
    %         end
    % %             errorbar(meanPRTA(:,cell),semPRTA(:,cell))
    %         line([2*Fspan,2*Fspan],[0,.5*nr]);
    %
    %         xlim([0 size(meanPRTA,1)])
    %         set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
    %         xlabel('Seconds')
    %         title(['PRTloc traces mouse ', MouseID]);
    %
    %         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTloc all traces ', num2str(tl*4),' seconds', '.jpg']);
    %         savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTloc all traces ', num2str(tl*4),' seconds','.fig']);
    %
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    %     splot=1;
    %     %plot all F with lines to show split
    %     for cell=fig_start:fig_end
    %         subplot(sqrt(num_square),sqrt(num_square),splot)
    %         if Fs==15.5/4
    %             t=0:1/Fs:(4*Fspan/Fs);%+1/Fs;
    %         else
    %             t=0:1/Fs:(4*Fspan/Fs);
    %         end
    t=0:(size(meanPRTA,1)-1);
    plot(t,meanPRTA(:),'b',t,(meanPRTA(:)+2*semPRTA(:)),'g:',t,(meanPRTA(:)-2*semPRTA(:)),'g:')
    hold on
    
    [pks,locs]=findpeaks(meanPRTA,'MinPeakDistance',Fs*3);
    [maxes,inds]=sort(pks,'descend');
    if ~isempty(inds)
    maxinfo{c}(1:4)=[pks(inds(1));locs(inds(1));pks(inds(min(2,length(inds))));locs(inds(min(2,length(inds))))];
    [pks,locs]=findpeaks(-meanPRTA(:),'MinPeakDistance',Fs*3);
    [maxes,inds]=sort(pks,'descend');
    maxinfo{c}(5:6)=[-pks(inds(1));locs(inds(1))];
    scatter(t(maxinfo{c}([2 4 6])),maxinfo{c}([1 3 5]));
    xes=maxinfo{c}([2 4 6]);
    maxinfo{c}([2 4 6])=(t(maxinfo{c}([2 4 6]))-ceil(size(meanPRTA,1)/2))/Fs;
    text(xes(1),maxinfo{c}(1)+.01,num2str(maxinfo{c}(2)))
    text(xes(2),maxinfo{c}(3)+.01,num2str(maxinfo{c}(4)))
    text(xes(3),maxinfo{c}(5)+.01,num2str(maxinfo{c}(6)))
    end
    set(gca,'XTick',[t(1:Fspan:end)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
    xlabel('Seconds')
    line([2*Fspan,2*Fspan],[min(meanPRTA(:)-3*semPRTA(:)),max(meanPRTA(:)+2*semPRTA(:))]);
    if c>1
    line([2*Fspan-avgrewdist,2*Fspan-avgrewdist],[min(meanPRTA(:)-3*semPRTA(:)),max(meanPRTA(:)+2*semPRTA(:))],'Color','r');
    line([2*Fspan+avgrewdist,2*Fspan+avgrewdist],[min(meanPRTA(:)-3*semPRTA(:)),max(meanPRTA(:)+2*semPRTA(:))],'Color','r');
    end
    xlim([0 size(PRTA,1)])
    
    title(['PRTA confidence mouse ',MouseID]);
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTloc confidence', num2str(tl*4),' seconds', env,'.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTloc confidence', num2str(tl*4),' seconds',env,'.fig']);
end
if exist('part','var')
    for c=2%:length(cPRTA)
            PRTA=cPRTA{c};
            tl=(size(PRTA,1)-1)/(4*Fs);
        
        import mlreportgen.dom.*;

        append(part,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'PRTloc heatmap', ' ', num2str(tl*4),' seconds',env,'.jpg'])));
        append(part,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'PRTloc confidence', num2str(tl*4),' seconds', env,'.jpg'])));

    end
end

%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         [~,prtainds]=max(meanPRTA);
%         [~,sortedinds]=sort(prtainds);
%         imagesc(meanPRTA(:,sortedinds)')
%         set(gca,'XTick',[0:Fspan:size(PRTA,1)],'XTickLabel',round([-size(PRTA,1)/2:Fspan:size(PRTA,1)/2]/Fs));
%         xlabel('Seconds')
%
%         saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap all', ' ', num2str(tl*4),' seconds', num2str(f), '.jpg']);
%         savefig([saveDir,'\',MouseID,'\', MouseID, 'PRTA heatmap all', ' ', num2str(tl*4),' seconds', num2str(f),'.fig']);
%
%
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
