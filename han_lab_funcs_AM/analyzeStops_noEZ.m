function excitedstop=analyzeStops_noEZ(mouse,varargin)
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
    Fs=4*ones(size(mouseNums));
    saveDir='F:\MA Data\Interneurons\PubFigs';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
percentstopped(length(mouseNums)*length(expnum))=0;
ii=1;
allstopped=[];
allgo=[];
excitedstop(6)=struct('Cells',cell(1));
for m=mouseNums
    for e=expnum{m}
        disp(['m',num2str(m),'e',num2str(e)])
        
        F1=mouse(m).Falls{e};
        y1=mouse(m).ybinned{e};
        f1=mouse(m).Forwards{e};
        r1=mouse(m).Rotations{e};
        y1=bsxfun(@minus,y1,min(y1,[],1));
        y1=bsxfun(@rdivide,y1,max(y1,[],1))*180;
        yind=(y1<170&y1>10);
        f2=nan(max(sum(yind,1)),size(F1,2));
        r2=f2;F2=f2;y2=f2;
        for c=1:size(F1,2)
            f2(:,c)=[f1(yind(:,c),c);nan(max(sum(yind,1))-sum(yind(:,c)),1)];
            r2(:,c)=[r1(yind(:,c),c);nan(max(sum(yind,1))-sum(yind(:,c)),1)];
            F2(:,c)=[F1(yind(:,c),c);nan(max(sum(yind,1))-sum(yind(:,c)),1)];
            y2(:,c)=[y1(yind(:,c),c);nan(max(sum(yind,1))-sum(yind(:,c)),1)];
        end
        
                %% Just forward motion
%                 [stFall,goFall,stopinds]=find_stops(F2,f2,Fs(m,e));
%         
%                 stoppedF=reshape(stFall,[],1);
%                 goingF=reshape(goFall,[],1);
%                 percentstopped(ii)=sum(isnan(goingF))/(length(goingF));
%                 allstopped=[stoppedF;allstopped];
%                 allgo=[goingF;allgo];
%         
%                 stopline=stopinds(:,1);
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 ygo=y2(:,1);
%                 ygo(stopline>0)=NaN;
%                 yst=y2(:,1);
%                 yst(stopline==0)=NaN;
%                 s1=subplot(2,1,1);
%                 plot(yst,'r')
%                 hold on
%                 plot(ygo,'g')
%                 s2=subplot(2,1,2);
%                 plot(f2(:,1).*(stopline==0),'g')
%                 hold on
%                 plot(f2(:,1).*(stopline>0),'r')
%                 linkaxes([s1,s2],'x')
%                 suptitle('Position and Speed of Stop and Go times w/o EZ ')
%         
%                 saveas(gca,[saveDir,'\', 'Stopped pos and speed wo EZ ',num2str(m),',',num2str(e), '.jpg']);
%                 savefig([saveDir,'\', 'Stopped pos and speed wo EZ ',num2str(m),',',num2str(e), '.fig']);
%         
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 mst=nanmean(stoppedF);
%                 mgo=nanmean(goingF);
%                 semst=1.96*nanstd(stoppedF)/sqrt(sum(isnan(goingF)));
%                 semgo=1.96*nanstd(goingF)/sqrt(sum(isnan(stoppedF)));
%                 subplot(1,2,1)
%                 hold on
%                 plot([mst;mgo],'k:','LineWidth',1,'Marker','x','MarkerSize',8,'MarkerFaceColor','k')
%                 errorbar([mst;mgo],[semst;semgo],'r.','LineWidth',.5)
%                 if ~isnan(nanmean(stoppedF))
%                 [p,h,stats]=ranksum(stoppedF,goingF,'tail','left');
%                 
%                 title(['Mean dF/F for Mouse ',num2str(m),' rec ', num2str(e), ' p= ',num2str(p)]);
%                 set(gca,'xtick',[1 2],'xticklabel', {'Stopped','Moving'})
%                 end
%                 subplot(1,2,2)
%                 histogram(stoppedF);
%                 hold on
%                 histogram(goingF);
%                 legend('Stopped','Moving')
%                 title(['Histogram of dF w/o EZ stopped percent=',num2str(percentstopped(ii))])
%         
%                 saveas(gca,[saveDir,'\', 'Stopped dF wo EZ ',num2str(m),',',num2str(e), '.jpg']);
%                 savefig([saveDir,'\', 'Stopped dF wo EZ ',num2str(m),',',num2str(e), '.fig']);
%         
%         
%                 figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%                 stoploc=find(stopline);
%                 usable=stoploc>(16*Fs(m,e))&stoploc<(length(F2)-(44*Fs(m,e)));
%                 %         stoploc=find(stopline)>(16*Fs(m,e))&find(stopline)<(44*Fs(m,e));
%                 if ~isempty(usable)
%                     usestop=stoploc(usable(1));
%                     if isempty(usestop)
%                         showstop=max(stoploc(1)-16*Fs(m,e),1):min(stoploc(1)+44*Fs(m,e),length(F2));
%                     else
%                         showstop=max(usestop-16*Fs(m,e),1):min(usestop+44*Fs(m,e),length(F2));
%                     end
%                     time=linspace(0,length(showstop)/Fs(m,e),length(showstop));
%                     for jj=1:(min(size(F2,2),8))
%                         subplot(4,2,jj)
%                         plot(time,F2(showstop,jj)+1)
%                         hold on
%                         plot(time,bsxfun(@rdivide,f2(showstop,jj),max(f2,[],1)))
%                         xlabel('Time (s)')
%                     end
%                     suptitle('Sample dF/F with speed')
%                     saveas(gca,[saveDir,'\', 'Sample dFF with speed wo EZ ',num2str(m),',',num2str(e), '.jpg']);
%                     savefig([saveDir,'\', 'Sample dFF with speed wo EZ ',num2str(m),',',num2str(e), '.fig']);
%                 end
        %% forward and rotation
        [stFall,goFall,stopinds]=find_stops(F2,sqrt(f2.^2+r2.^2),Fs(m,e));
        
        stoppedF=reshape(stFall,[],1);
        goingF=reshape(goFall,[],1);
        percentstopped(ii)=sum(isnan(goingF))/(length(goingF));
        allstopped=[stoppedF;allstopped];
        allgo=[goingF;allgo];
        stopline=stopinds(:,1);
        disp(['Mouse: ', num2str(m),' Rec: ', num2str(e),' stop>go ', num2str(find(nanmean(stFall)>nanmean(goFall)))])
        excitedstop(m).Cells{e}=find(nanmean(stFall)>nanmean(goFall));
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        ygo=y2(:,1);
        ygo(stopline>0)=NaN;
        yst=y2(:,1);
        yst(stopline==0)=NaN;
        s1=subplot(2,1,1);
        plot(yst,'r')
        hold on
        plot(ygo,'g')
        s2=subplot(2,1,2);
        plot(sqrt(f2(:,1).^2+r2(:,1).^2).*(stopline==0),'g')
        hold on
        plot(sqrt(f2(:,1).^2+r2(:,1).^2).*(stopline>0),'r')
        linkaxes([s1,s2],'x')
        suptitle('Position and Speed of Stop and Go times with Rotation w/o EZ ')
        if ~exist([saveDir,'\', 'Stopped Pos and Speed\'],'dir')
            mkdir([saveDir,'\', 'Stopped Pos and Speed\'])
        end
        saveas(gca,[saveDir,'\', 'Stopped Pos and Speed\','Stopped pos and speed with Rotation  wo EZ  ',num2str(m),',',num2str(e), '.jpg']);
        savefig([saveDir,'\',  'Stopped Pos and Speed\','Stopped pos and speed with Rotation wo EZ ',num2str(m),',',num2str(e), '.fig']);
        
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        mst=nanmean(stoppedF);
        mgo=nanmean(goingF);
        semst=1.96*nanstd(stoppedF)/sqrt(sum(isnan(goingF)));
        semgo=1.96*nanstd(goingF)/sqrt(sum(isnan(stoppedF)));
        subplot(1,2,1)
        hold on
        errorbar([mst;mgo],[semst;semgo],'r.','LineWidth',.5)
        plot([mst;mgo],'w:','LineWidth',1,'Marker','o','MarkerSize',8,'MarkerFaceColor','k')
        testF=stoppedF;
        testF(isnan(stoppedF))=[];
        if ~isempty(testF)
            [p,h,stats]=ranksum(stoppedF,goingF);
        end
        if exist('p','var')
        title(['Mean dF/F for Mouse ',num2str(m),' rec ', num2str(e), ' p= ',num2str(p)]);
       end
        set(gca,'xtick',[1 2],'xticklabel', {'Stopped','Moving'})
        
        subplot(1,2,2)
        hold on
        histogram(goingF,'FaceAlpha',.2,'EdgeColor','none');
        
        histogram(stoppedF);
        
        legend('Moving','Stopped')
        title(['Histogram of dF with rotation w/o EZ stopped percent=',num2str(percentstopped(ii))])
        
        if ~exist([saveDir,'\', 'Stopped df with Rotation\'],'dir')
            mkdir([saveDir,'\', 'Stopped df with Rotation\'])
        end
        
        saveas(gca,[saveDir,'\', 'Stopped df with Rotation\', 'Stopped dF with Rotation  wo EZ ',num2str(m),',',num2str(e), '.jpg']);
        savefig([saveDir,'\', 'Stopped df with Rotation\', 'Stopped dF with Rotation wo EZ ',num2str(m),',',num2str(e), '.fig']);
        
        if ~exist([saveDir,'\', 'Sample dFF with speed\'],'dir')
            mkdir([saveDir,'\', 'Sample dFF with speed\'])
        end
        
        stoploc=find(stopline);
        usable=stoploc>(16*Fs(m,e))&stoploc<(length(F2)-(44*Fs(m,e)));
        if ~isempty(usable)
            usestop=stoploc(usable(1));
            if isempty(usestop)
                showstop=max(stoploc(1)-16*Fs(m,e),1):min(stoploc(1)+44*Fs(m,e),length(F2));
            else
                showstop=max(usestop-16*Fs(m,e),1):min(usestop+44*Fs(m,e),length(F2));
            end
            time=linspace(0,length(F2)/Fs(m,e),length(F2));
            for f=1:ceil(size(F2,2)/8)
                figure('units','normalized', 'Position', [.01 .05 .98 .87]);

                for jj=(1+(f-1)*8):(min(size(F2,2),f*8))
                    subplot(4,2,jj-(f-1)*8)
                    plot(time,F2(:,jj)+1)
                    hold on
                    plot(time,sqrt(f2(:,jj).^2+r2(:,jj).^2)/max(sqrt(f2(:,jj).^2+r2(:,jj).^2)))
                    xlabel('Time (s)')
                    xlim([showstop(1)/Fs(m,e) showstop(end)/Fs(m,e)])
                    title(['Cell: ', num2str(jj)])
                end
                suptitle('Sample dF/F with speed')
                
                
                saveas(gca,[saveDir,'\', 'Sample dFF with speed\','Sample dFF with speed rotation wo EZ ',num2str(m),',',num2str(e), 'Fig',num2str(f),'.jpg']);
                savefig([saveDir,'\', 'Sample dFF with speed\','Sample dFF with speed rotation wo EZ ',num2str(m),',',num2str(e), 'Fig',num2str(f), '.fig']);
                
            end
        else
            show=(50*Fs(m,e)):(110*Fs(m,e));
            time=linspace(0,length(F2)/Fs(m,e),length(F2));
             for f=1:ceil(size(F2,2)/8)
                figure('units','normalized', 'Position', [.01 .05 .98 .87]);

                for jj=(1+(f-1)*8):(min(size(F2,2),f*8))
                    subplot(4,2,jj-(f-1)*8)
                    plot(time,F2(:,jj)+1)
                    hold on
                    plot(time,sqrt(f2(:,jj).^2+r2(:,jj).^2)/max(sqrt(f2(:,jj).^2+r2(:,jj).^2)))
                    xlabel('Time (s)')
                    xlim([show(1)/Fs(m,e) show(end)/Fs(m,e)])
                    title(['Cell: ', num2str(jj)])
                end
                suptitle('Sample dF/F with speed')
                
                
                saveas(gca,[saveDir,'\', 'Sample dFF with speed\','Sample dFF with speed rotation wo EZ ',num2str(m),',',num2str(e), 'Fig',num2str(f),'.jpg']);
                savefig([saveDir,'\', 'Sample dFF with speed\','Sample dFF with speed rotation wo EZ ',num2str(m),',',num2str(e), 'Fig',num2str(f), '.fig']);
                
            end
            
            
            
        end
        
        close all
    end
end

figure('units','normalized', 'Position', [.01 .05 .98 .87]);
mst=nanmean(allstopped);
mgo=nanmean(allgo);
semst=1.96*nanstd(allstopped)/sqrt(sum(isnan(allgo)));
semgo=1.96*nanstd(allgo)/sqrt(sum(isnan(allstopped)));
subplot(1,2,1)
hold on
errorbar([mst;mgo],[semst;semgo],'r.','LineWidth',.5)
plot([mst;mgo],'w:','LineWidth',1,'Marker','o','MarkerSize',8,'MarkerFaceColor','k')
[p,h,stats]=ranksum(allstopped,allgo);
title(['Mean dF/F for all mice', ' p= ',num2str(p)]);
set(gca,'xtick',[1 2],'xticklabel', {'Stopped','Moving'})
subplot(1,2,2)
hold on       
histogram(allgo,'FaceAlpha',.2,'EdgeColor','none');
histogram(allstopped);
legend('Stopped','Moving')
title('Histogram of dF')

saveas(gca,[saveDir,'\', 'Stopped dF',num2str(m),',',num2str(e), '.jpg']);
savefig([saveDir,'\', 'Stopped dF',num2str(m),',',num2str(e), '.fig']);

