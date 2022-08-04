function analyzeStops(mouse,varargin)
if nargin>1
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>4
        saveDir=varargin{4};
    else     saveDir='F:\MA Data\Interneurons\PubFigs';
    end
    
else
    mouseNums=1:length(mouse);
    for m=mouseNums
        expnum{m}=1:length(mouse(m).Falls);
    end
    Fs=4*ones(size(mouseNums),max(expnum{1}));
    saveDir='F:\MA Data\Interneurons\PubFigs';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
percentstopped(length(mouseNums)*length(expnum))=0;
ii=1;
allstopped=[];
allgo=[];
for m=mouseNums
    for e=expnum{m}
        disp(['m',num2str(m),'e',num2str(e)])
        Fall=mouse(m).Falls{e};
        ybinned=mouse(m).ybinned{e};
        forward=mouse(m).Forwards{e};
        rotation=mouse(m).Rotations{e};
        rewards=mouse(m).rewards{e};
        rewratio=(length(rewards)/size(Fall,1));
        rewsmall2=zeros(length(Fall),1);
        
        num_cells=size(Fall,2);
        if sum(rewards)>0
        for jj=1:length(Fall)
            rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            [~,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs(m,e))*2);
        end
        else
            rewsmall2=zeros(size(Fall,1),1);
        end
        rewsmall=zeros(size(rewsmall2));
        rewsmall((rewlocs))=1;
        
        %% Just forward motion
%         [stFall,goFall,stopinds]=find_stops(Fall,forward,Fs(m,e));
%         stopline=logical(stopinds(:,1));
%         
%         stoppedF=reshape(stFall,[],1);
%         goingF=reshape(goFall,[],1);
%         percentstopped(ii)=sum(isnan(goingF))/(length(goingF));
%         allstopped=[stoppedF;allstopped];
%         allgo=[goingF;allgo];
%         
%         
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         ygo=ybinned(:,1);
%         ygo(stopline)=NaN;
%         yst=ybinned(:,1);
%         yst(~stopline)=NaN;
%         s1=subplot(2,1,1);
%         plot(yst,'r')
%         hold on
%         plot(ygo,'g')
%         s2=subplot(2,1,2);
%         plot(forward(:,1).*(~stopline),'g')
%         hold on
%         plot(forward(:,1).*(stopline),'r')
%         linkaxes([s1,s2],'x')
%         suptitle('Position and Speed of Stop and Go times')
%         
%         saveas(gca,[saveDir,'\', 'Stopped pos and speed',num2str(m),',',num2str(e), '.jpg']);
%         savefig([saveDir,'\', 'Stopped pos and speed',num2str(m),',',num2str(e), '.fig']);
%         
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         mst=nanmean(stoppedF);
%         mgo=nanmean(goingF);
%         semst=1.96*nanstd(stoppedF)/sqrt(sum(isnan(goingF)));
%         semgo=1.96*nanstd(goingF)/sqrt(sum(isnan(stoppedF)));
%         subplot(1,2,1)
%         hold on
%         plot([mst;mgo],'k:','LineWidth',1,'Marker','x','MarkerSize',8,'MarkerFaceColor','k')
%         errorbar([mst;mgo],[semst;semgo],'r.','LineWidth',.5)
%         [h,p,ci,stats]=ttest2(stoppedF,goingF,'varType','unequal');
%         title(['Mean dF/F for Mouse ',num2str(m),' rec ', num2str(e), ' p= ',num2str(p)]);
%         set(gca,'xtick',[1 2],'xticklabel', {'Stopped','Moving'})
%         subplot(1,2,2)
%         histogram(stoppedF);
%         hold on
%         histogram(goingF);
%         legend('Stopped','Moving')
%         title(['Histogram of dF percent=',num2str(percentstopped(ii))])
%         
%         saveas(gca,[saveDir,'\', 'Stopped dF',num2str(m),',',num2str(e), '.jpg']);
%         savefig([saveDir,'\', 'Stopped dF',num2str(m),',',num2str(e), '.fig']);
%         
%         
%         
%         stoploc=find(stopline);
%         usable=stoploc>(16*Fs(m,e))&stoploc<(length(Fall)-(44*Fs(m,e)));
%         %         stoploc=find(stopline)>(16*Fs(m,e))&find(stopline)<(44*Fs(m,e));
%         if ~isempty(usable)
% %             usestop=stoploc(usable(1));
%             usestop=60*4;
%             if isempty(usestop)
%                 showstop=max(stoploc(1)-16*Fs(m,e),1):min(stoploc(1)+44*Fs(m,e),length(Fall));
%             else
%                 showstop=max(usestop-16*Fs(m,e),1):min(usestop+64*Fs(m,e),length(Fall));
%             end
%             time=linspace(0,length(Fall)/Fs(m,e),length(Fall));
%             for f=1:ceil(size(Fall,2)/8)
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 
%                 for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
%                     s(jj)=subplot(4,2,jj-(f-1)*8);
%                     plot(time,Fall(:,jj)+1)
%                     hold on
%                     plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
%                     plot(time,2*(rewsmall2));
%                     xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
%                     xlabel('Time (s)')
%                 end
%                 suptitle('Sample dF/F with speed')
%                 linkaxes(s,'x')
% 
%                 saveas(gca,[saveDir,'\', 'Sample dFF with speed  ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
%                 savefig([saveDir,'\', 'Sample dFF with speed  ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);
%                
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 
%                 for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
%                     s(jj)=subplot(4,2,jj-(f-1)*8);
%                     plot(time,Fall(:,jj)+1)
%                     hold on
%                     plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
% %                     plot(time,2*rewsmall(:));
%                     xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
%                     xlabel('Time (s)')
%                 end
%                 linkaxes(s,'x')
% 
%                 suptitle('Sample dF/F with speed')
%                 saveas(gca,[saveDir,'\', 'Sample dFF with speed No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
%                 savefig([saveDir,'\', 'Sample dFF with speed No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);
% 
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 
%                 for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
%                     plot(time,Fall(:,jj)+1)
%                     hold on
%                     plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
%                     plot(time,2*(rewsmall2));
%                     xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
%                     xlabel('Time (s)')
%                 end
%                 suptitle('Sample dF/F with speed')
%                 saveas(gca,[saveDir,'\', 'Sample dFF with speed Overlay ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
%                 savefig([saveDir,'\', 'Sample dFF with speed Overlay  ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);
%             
%             end
%         else
%             show=(60*Fs(m,e)):(120*Fs(m,e));
%             time=linspace(0,length(Fall)/Fs(m,e),length(Fall));
%             for f=1:ceil(size(Fall,2)/8)
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 
%                 for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
%                     s(jj)=subplot(4,2,jj-(f-1)*8);
%                     plot(time,Fall(:,jj)+1)
%                     hold on
%                     plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
%                     plot(time,2*rewsmall(:));
%                     xlim([show(1)/Fs(m,e) show(end)/Fs(m,e)])
%                     xlabel('Time (s)')
%                 end
%                 linkaxes(s,'x')
% 
%                 suptitle('Sample dF/F with speed')
%                 saveas(gca,[saveDir,'\', 'Sample dFF with speed  ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
%                 savefig([saveDir,'\', 'Sample dFF with speed  ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);
% 
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 
%                 for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
%                     s(jj)=subplot(4,2,jj-(f-1)*8);
%                     plot(time,Fall(:,jj)+1)
%                     hold on
%                     plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
% %                     plot(time,2*rewsmall(:));
%                     xlim([show(1)/Fs(m,e) show(end)/Fs(m,e)])
%                     xlabel('Time (s)')
%                 end
%                 linkaxes(s,'x')
% 
%                 suptitle('Sample dF/F with speed')
%                 saveas(gca,[saveDir,'\', 'Sample dFF with speed No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
%                 savefig([saveDir,'\', 'Sample dFF with speed No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);
% 
% 
% 
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 
%                 for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
%                     plot(time,Fall(:,jj)+1)
%                     hold on
%                     plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
%                     plot(time,2*(rewsmall2));
%                     xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
%                     xlabel('Time (s)')
%                 end
%                 suptitle('Sample dF/F with speed')
%                 saveas(gca,[saveDir,'\', 'Sample dFF with speed Overlay ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
%                 savefig([saveDir,'\', 'Sample dFF with speed Overlay  ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);
%             end
%         end
        
        
        %forward and rotation
        total=sqrt(forward.^2+rotation.^2);
        [stFall,goFall,stopinds]=find_stops(Fall,total,Fs(m,e));
        
        stoppedF=reshape(stFall,[],1);
        goingF=reshape(goFall,[],1);
        percentstopped(ii)=sum(isnan(goingF))/(length(goingF));
        allstopped=[stoppedF;allstopped];
        allgo=[goingF;allgo];
        
        stopline=logical(stopinds(:,1));
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        ygo=ybinned(:,1);
        ygo(stopline)=NaN;
        yst=ybinned(:,1);
        yst(~stopline)=NaN;
        s1=subplot(2,1,1);
        plot(yst,'r')
        hold on
        plot(ygo,'g')
        s2=subplot(2,1,2);
        plot(total(:,1).*(~stopline),'g')
        hold on
        plot(total(:,1).*(stopline),'r')
        linkaxes([s1,s2],'x')
        suptitle('Position and Speed of Stop and Go times with Rotation')
        
        saveas(gca,[saveDir,'\', 'Stopped pos and speed with Rotation ',num2str(m),',',num2str(e), '.jpg']);
        savefig([saveDir,'\', 'Stopped pos and speed with Rotation ',num2str(m),',',num2str(e), '.fig']);
        
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        mst=nanmean(stoppedF);
        mgo=nanmean(goingF);
        semst=1.96*nanstd(stoppedF)/sqrt(sum(isnan(goingF)));
        semgo=1.96*nanstd(goingF)/sqrt(sum(isnan(stoppedF)));
        subplot(1,2,1)
        hold on
        plot([mst;mgo],'k:','LineWidth',1,'Marker','x','MarkerSize',8,'MarkerFaceColor','k')
        errorbar([mst;mgo],[semst;semgo],'r.','LineWidth',.5)
        [h,p,ci,stats]=ttest2(stoppedF,goingF,'varType','unequal');
        title(['Mean dF/F for Mouse ',num2str(m),' rec ', num2str(e), ' p= ',num2str(p)]);
        set(gca,'xtick',[1 2],'xticklabel', {'Stopped','Moving'})
        subplot(1,2,2)
        histogram(stoppedF);
        hold on
        histogram(goingF);
        legend('Stopped','Moving')
        title(['Histogram of dF with Rotation percent=',num2str(percentstopped(ii))])
        
        saveas(gca,[saveDir,'\', 'Stopped dF with Rotation',num2str(m),',',num2str(e), '.jpg']);
        savefig([saveDir,'\', 'Stopped dF with Rotation',num2str(m),',',num2str(e), '.fig']);
        
        stoploc=find(stopline);
        usable=stoploc>(16*Fs(m,e))&stoploc<(length(Fall)-(44*Fs(m,e)));
        if ~isempty(usable)
            usestop=stoploc(usable(1));
%             usestop=1270*4;
            if isempty(usestop)
                showstop=max(stoploc(1)-16*Fs(m,e),1):min(stoploc(1)+44*Fs(m,e),length(Fall));
            else
                showstop=max(usestop-16*Fs(m,e),1):min(usestop+64*Fs(m,e),length(Fall));
            end
            time=linspace(0,length(Fall)/Fs(m,e),length(Fall));
            for f=1:ceil(size(Fall,2)/8)
                figure('units','normalized', 'Position', [.01 .05 .98 .87]);
                1
                for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
                    s(jj)=subplot(4,2,jj-(f-1)*8);
                    plot(time,Fall(:,jj)+1)
                    hold on
                    plot(time,bsxfun(@rdivide,total(:,jj),max(total(:),[],1)))
                    plot(time,2*rewsmall(:));
                    xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
                    
                    xlabel('Time (s)')
                end
                linkaxes(s,'x')

                suptitle('Sample dF/F with speed')
                saveas(gca,[saveDir,'\', 'Sample dFF with speed rotation ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
                savefig([saveDir,'\', 'Sample dFF with speed rotation  ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);

                figure('units','normalized', 'Position', [.01 .05 .98 .87]);
                2
                for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
                    s(jj)=subplot(4,2,jj-(f-1)*8);
                    plot(time,Fall(:,jj)+1)
                    hold on
                    plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
%                     plot(time,2*rewsmall(:));
                    xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
                    xlabel('Time (s)')
                end
                linkaxes(s,'x')

                suptitle('Sample dF/F with speed')
                saveas(gca,[saveDir,'\', 'Sample dFF with speed rotation No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
                savefig([saveDir,'\', 'Sample dFF with speed rotation No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);

                figure('units','normalized', 'Position', [.01 .05 .98 .87]);
                3
                for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
                    plot(time,Fall(:,jj)+1.5)
                    hold on
%                     plot(time,2*(rewsmall));
                    xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
                    xlabel('Time (s)')
                end
                    plot(time,bsxfun(@rdivide,forward(:,1),max(forward(:),[],1)))

                suptitle('Sample dF/F with speed')
                saveas(gca,[saveDir,'\', 'Sample dFF with speed rotation Overlay No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
                savefig([saveDir,'\', 'Sample dFF with speed rotation Overlay No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);

            end
            
        else
            show=(60*Fs(m,e)):(120*Fs(m,e));
            time=linspace(0,length(Fall)/Fs(m,e),length(Fall));
            
            for f=1:ceil(size(Fall,2)/8)
                figure('units','normalized', 'Position', [.01 .05 .98 .87]);
                
                for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
                    s(jj)=subplot(4,2,jj-(f-1)*8);
                    plot(time,Fall(:,jj)+1)
                    hold on
                    plot(time,bsxfun(@rdivide,total(:,jj),max(total(:),[],1)))
                    plot(time,2*rewsmall(:));
                    xlim([show(1)/Fs(m,e) show(end)/Fs(m,e)])
                    xlabel('Time (s)')
                end
                linkaxes(s,'x')
                suptitle('Sample dF/F with speed')
                saveas(gca,[saveDir,'\', 'Sample dFF with speed rotation ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
                savefig([saveDir,'\', 'Sample dFF with speed rotation ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);


                figure('units','normalized', 'Position', [.01 .05 .98 .87]);
                
                for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
                    s(jj)=subplot(4,2,jj-(f-1)*8);
                    plot(time,Fall(:,jj)+1)
                    hold on
                    plot(time,bsxfun(@rdivide,forward(:,jj),max(forward(:),[],1)))
%                     plot(time,2*rewsmall(:));
                    xlim([show(1)/Fs(m,e) show(end)/Fs(m,e)])
                    xlabel('Time (s)')
                end
                linkaxes(s,'x')

                suptitle('Sample dF/F with speed')
                saveas(gca,[saveDir,'\', 'Sample dFF with speed rotation No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
                savefig([saveDir,'\', 'Sample dFF with speed rotation No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.fig']);

                figure('units','normalized', 'Position', [.01 .05 .98 .87]);
                
                for jj=(1+(f-1)*8):(min(size(Fall,2),f*8))
                    plot(time,Fall(:,jj)+1.5)
                    hold on
%                     plot(time,2*(rewsmall));
                    xlim([show(1)/Fs(m,e) show(end)/Fs(m,e)])
                    xlabel('Time (s)')
                end
                    plot(time,bsxfun(@rdivide,forward(:,1),max(forward(:),[],1)))

                suptitle('Sample dF/F with speed')
                saveas(gca,[saveDir,'\', 'Sample dFF with speed rotation Overlay No rew ',num2str(m),',',num2str(e),'Fig',num2str(f), '.jpg']);
                savefig
            end
        end
    end
end


figure('units','normalized', 'Position', [.01 .05 .98 .87]);
mst=nanmean(allstopped);
mgo=nanmean(allgo);
semst=1.96*nanstd(allstopped)/sqrt(sum(isnan(allgo)));
semgo=1.96*nanstd(allgo)/sqrt(sum(isnan(allstopped)));
subplot(1,2,1)
hold on
plot([mst;mgo],'k:','LineWidth',1,'Marker','x','MarkerSize',8,'MarkerFaceColor','k')
errorbar([mst;mgo],[semst;semgo],'r.','LineWidth',.5)
[h,p,ci,stats]=ttest2(allstopped,allgo,'varType','unequal');
title(['Mean dF/F for all mice', ' p= ',num2str(p)]);
set(gca,'xtick',[1 2],'xticklabel', {'Stopped','Moving'})
subplot(1,2,2)
histogram(allstopped);
hold on
histogram(allgo);
legend('Stopped','Moving')
title('Histogram of dF')

saveas(gca,[saveDir,'\', 'Stopped dF',num2str(m),',',num2str(e), '.jpg']);
savefig([saveDir,'\', 'Stopped dF',num2str(m),',',num2str(e), '.fig']);