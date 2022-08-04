function [collect_means,startinfo,stopinfo]=xcorr_speed_change(mouse,varargin)
if nargin>1
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>4
        saveDir=varargin{4};
    else     saveDir='F:\MA Data\Interneurons\PubFigs\Xcorr\';
    end
    
else
    mouseNums=1:length(mouse);
    for m=mouseNums
        expnum{m}=1:length(mouse(m).Falls);
    end
    Fs=4*ones(size(mouseNums));
    saveDir='F:\MA Data\Interneurons\PubFigs\Xcorr\';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
collect_means(length(mouse))=struct('starting',cell(1),'stopping',cell(1));
parfor m=mouseNums
    for e=expnum{m}
        allstarting=[];
        allstopping=[];
        allstartingxc=[];
        allstoppingxc=[];
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
        
        binsize=5*Fs(m,e);
        [~,~,stopinds]=find_stops(F2,sqrt(f2.^2+r2.^2),Fs(m,e));
        dstopinds=[zeros(1,size(stopinds,2));diff(stopinds)];
        beginstop=find(dstopinds==1);
        [bsr,bsc]=ind2sub(size(dstopinds),beginstop);
        stoppingF=zeros([binsize+1+Fs(m,e),ceil(max(size(bsr,1))/max(bsc)),max(bsc)]);
        stoppingV=stoppingF;
        curcol=1;
        colcount=1;
        for ii=1:length(bsr)
            if bsc(ii)~=curcol
                curcol=bsc(ii);
                colcount=1;
            end
            start=max((bsr(ii)-binsize),1);
            stop=min(bsr(ii)+Fs(m,e),length(F2));
            tempF=F2(start:stop,bsc(ii));
            tempV=f2(start:stop,bsc(ii));
            if (bsr(ii)-binsize)<1
                tempF=[nan(size(stoppingF,1)-size(tempF,1),1); tempF];
                tempV=[nan(size(stoppingF,1)-size(tempV,1),1); tempV];
            elseif (bsr(ii)+Fs(m,e))>length(F2)
                tempF=[tempF;nan(size(stoppingF,1)-size(tempF,1),1)];
                tempV=[tempV;an(size(stoppingF,1)-size(tempV,1),1)];
            end
            stoppingF(:,colcount,bsc(ii))=tempF;
            stoppingV(:,colcount,bsc(ii))=tempV;
            
            colcount=colcount+1;
        end
        
        endstop=find(dstopinds==-1);
        [esr,esc]=ind2sub(size(dstopinds),endstop);
        startingF=zeros([binsize+1+Fs(m,e),ceil(length(esr)/max(esc)),max(esc)]);
        startingV=startingF;
        curcol=1;
        colcount=1;
        for ii=1:length(esr)
            if esc(ii)~=curcol
                curcol=esc(ii);
                colcount=1;
            end
            start=max(esr(ii)-Fs(m,e),1);
            stop=min(esr(ii)+binsize,length(F2));
            tempF=F2(start:stop,esc(ii));
            tempV=f2(start:stop,esc(ii));
            if (esr(ii)+binsize)>length(F2)
                tempF=[tempF;nan(size(startingF,1)-size(tempF,1),1)];
                tempV=[tempV;nan(size(startingF,1)-size(tempV,1),1)];
            elseif (esr(ii)-Fs(m,e))<1
                tempF=[nan(size(stoppingF,1)-size(tempF,1),1); tempF];
                tempV=[nan(size(stoppingF,1)-size(tempV,1),1); tempV];
            end
            startingF(:,colcount,esc(ii))=tempF;
            startingV(:,colcount,esc(ii))=tempV;
            colcount=colcount+1;
            
        end
        timestarting=linspace(-1,5,binsize+1+Fs(m,e));
        timestopping=linspace(-5,1,binsize+1+Fs(m,e));
        collect_means(m).starting{e}=[];
        collect_means(m).stopping{e}=[];
        startinfo=zeros(2,size(F2,2));
        stopinfo=zeros(2,size(F2,2));
        if ~isempty(esr)
            for cellnum=1:size(F2,2)
                if size(startingF,3)>=cellnum
                    mstartingF=nanmean(squeeze(startingF(:,:,cellnum)),2);
                    semstartingF=nanstd(squeeze(startingF(:,:,cellnum)),0,2)/sqrt(size(startingF,2));
                    mstartingV=nanmean(squeeze(startingV(:,:,cellnum)),2);
                    semstartingV=nanstd(squeeze(startingV(:,:,cellnum)),0,2)/sqrt(size(startingF,2));
                    [xcstart,lagstart]=xcorr(mstartingF-mean(mstartingF),mstartingV-mean(mstartingV),'coeff');
                    mstoppingF=nanmean(squeeze(stoppingF(:,:,cellnum)),2);
                    semstoppingF=nanstd(squeeze(stoppingF(:,:,cellnum)),0,2)/sqrt(size(startingF,2));
                    mstoppingV=nanmean(squeeze(stoppingV(:,:,cellnum)),2);
                    semstoppingV=nanstd(squeeze(stoppingV(:,:,cellnum)),0,2)/sqrt(size(startingF,2));
                    [xcstop,lagstop]=xcorr(mstoppingF-mean(mstoppingF),mstoppingV-mean(mstoppingV),'coeff');
                    
                    collect_means(m).starting{e}=[collect_means(m).starting{e},mstartingF];
                    collect_means(m).stopping{e}=[collect_means(m).stopping{e},mstoppingF];
                    
                    figure('units','normalized', 'Position', [.01 .05 .49 .87]);
                    s1=subplot(3,2,1);
                    plot(timestarting,mstartingF)
                    hold on
                    plot(timestarting,mstartingF+1.96*semstartingF,'g--')
                    plot(timestarting,mstartingF-1.96*semstartingF,'g--')
                    line([0,0],[min(mstartingF)-3*min(semstartingF),max(mstartingF)+3*max(semstartingF)])
                    if ~isnan(min(mstartingF)-3*min(semstartingF))
                        ylim([min(mstartingF)-3*min(semstartingF),max(max(mstartingF)+3*max(semstartingF),.1)])
                    end
                    title(['Activity from Stop, Cell ',num2str(cellnum)])
                    s2=subplot(3,2,3);
                    plot(timestarting,mstartingV)
                    hold on
                    plot(timestarting,mstartingV+1.96*semstartingV,'g--')
                    plot(timestarting,mstartingV-1.96*semstartingV,'g--')
                    line([0,0],[min(mstartingV)-3*min(semstartingV),max(mstartingV)+3*max(semstartingV)])
                    title(['Speed from Stop, Cell ',num2str(cellnum)])
                    s3=subplot(3,2,2);
                    plot(timestopping,mstoppingF)
                    hold on
                    plot(timestopping,mstoppingF+1.96*semstoppingF,'g--')
                    plot(timestopping,mstoppingF-1.96*semstoppingF,'g--')
                    line([0,0],[min(mstoppingF)-3*min(semstoppingF),max(mstoppingF)+3*max(semstoppingF)])
                    if ~isnan(min(mstoppingF)-3*min(semstoppingF))
                        ylim([min(mstoppingF)-3*min(semstoppingF),max(max(mstoppingF)+3*max(semstoppingF),.1)])
                    end
                    title(['Activity towards stop, Cell ',num2str(cellnum)])
                    s4=subplot(3,2,4);
                    plot(timestopping,mstoppingV)
                    hold on
                    plot(timestopping,mstoppingV+1.96*semstoppingV,'g--')
                    plot(timestopping,mstoppingV-1.96*semstoppingV,'g--')
                    line([0,0],[min(mstoppingV)-3*min(semstoppingV),max(max(mstoppingV)+3*max(semstoppingV),.1)])
                    title(['Speed towards stop, Cell ',num2str(cellnum)])
                    
                    subplot(3,2,5)
                    plot(lagstart/Fs(m,e),xcstart)
                    hold on
                    line([0 0], [min(xcstart),max(xcstart)])
                    xlabel('Seconds')
                    [startinfo(1,cellnum),startinfo(2,cellnum)]=max(xcstart);
                    startinfo(2,cellnum)=lagstart(startinfo(2,cellnum))/Fs(m,e);
                    scatter(startinfo(2,cellnum),startinfo(1,cellnum),'k+')
                    text(startinfo(2,cellnum),startinfo(1,cellnum)+.02*startinfo(1,cellnum),num2str(startinfo(2,cellnum)))
                    
                    subplot(3,2,6)
                    plot(lagstop/Fs(m,e),xcstop)
                    hold on
                    line([0 0], [min(xcstop),max(xcstop)])
                    xlabel('Seconds')
                    [stopinfo(1,cellnum),stopinfo(2,cellnum)]=max(xcstop);
                    stopinfo(2,cellnum)=lagstop(stopinfo(2,cellnum))/Fs(m,e);
                    scatter(stopinfo(2,cellnum),stopinfo(1,cellnum),'k+')
                    text(stopinfo(2,cellnum),stopinfo(1,cellnum)+.02*stopinfo(1,cellnum),num2str(stopinfo(2,cellnum)))
                    
                    linkaxes([s1,s2],'x')
                    linkaxes([s3,s4],'x')
                    
                    saveas(gca,[saveDir,'\', 'XCorr from stop',num2str(m),',',num2str(e),' Cell ',num2str(cellnum), '.jpg']);
                    savefig([saveDir,'\', 'XCorr from stop',num2str(m),',',num2str(e), ' Cell ',num2str(cellnum), '.fig']);
                    disp(['Mouse ',num2str(m),' Rec ',num2str(e)]);
                    allstarting=[allstarting,(mstartingF)];
                    allstartingxc=[allstartingxc,(xcstart)];
                    allstopping=[allstopping,(mstoppingF)];
                    allstoppingxc=[allstoppingxc,(xcstop)];
                    
                else
                    collect_means(m).starting{e}=[collect_means(m).starting{e},nan(size(startingF,1),1)];
                    collect_means(m).stopping{e}=[collect_means(m).stopping{e},nan(size(startingF,1),1)];
                end
            end
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            subplot(2,1,1)
            imagesc(allstarting')
            title('All Cell activity from stop')
            set(gca,'XTick',1:Fs(m,e):size(allstarting,1),'XTickLabel',{'-1','0','1','2','3','4','5'})
            subplot(2,1,2)
            imagesc(allstopping')
            set(gca,'XTick',1:Fs(m,e):size(allstarting,1),'XTickLabel',{'5','4','3','2','1','0','-1'})
            title('All Cell activity towards stop')
            saveas(gca,[saveDir,'\', 'XCorr from stop',num2str(m),',',num2str(e),' All Cells','.jpg']);
            savefig([saveDir,'\', 'XCorr from stop',num2str(m),',',num2str(e), ' All Cells ','.fig']);
            
            close all
        else
            collect_means(m).starting{e}=[collect_means(m).starting{e},nan(size(startingF,1),size(F2,2))];
            collect_means(m).stopping{e}=[collect_means(m).stopping{e},nan(size(startingF,1),size(F2,2))];
        end
    end
end