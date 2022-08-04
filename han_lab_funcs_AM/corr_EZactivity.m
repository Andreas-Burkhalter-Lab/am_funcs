function collect_means=corr_EZactivity(mouse,varargin)
if nargin>1
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>4
        saveDir=[varargin{4},'\EZactivity\'];
    else     saveDir='F:\MA Data\Interneurons\PubFigs\EZactivity\';
    end
    if nargin>5
        env=varargin{5};
    else env=' ';
    end
    
else
    mouseNums=1:length(mouse);
    for m=mouseNums
        expnum{m}=1:length(mouse(m).Falls);
    end
    Fs=4*ones(size(mouseNums));
    saveDir='F:\MA Data\Interneurons\PubFigs\EZactivity\';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
collect_means(length(mouse))=struct('starting',cell(1),'stopping',cell(1));
parfor m=mouseNums
    for e=expnum{m}
        m
        e
        
        allupwardsF=[];
        alldownwardsF=[];
        allupwardsV=[];
        alldownwardsV=[];
        F1=mouse(m).Falls{e};
        y1=mouse(m).ybinned{e};
        f1=mouse(m).Forwards{e};
        r1=mouse(m).Rotations{e};
        y1=bsxfun(@minus,y1,min(y1,[],1));
        y1=bsxfun(@rdivide,y1,max(y1,[],1))*180;
        for i=1:size(y1,2)
            y1(:,i)=smooth(y1(:,i));
        end
        rewards=mouse(m).rewards{e};
        rewratio=(length(rewards)/size(F1,1));
        rewsmall2=zeros(length(F1),1);
        num_cells=size(F1,2);
        for jj=1:length(F1)
            rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
        end
        [~,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs(m,e))*2);
        
        downcrosses=find((diff(y1(:,1)-10<0))==1);
        downleaving=find((diff(y1(:,1)-10<0))==-1);
        rewardsdown=downcrosses(ismembertol(downcrosses,rewlocs,round(Fs(m,e))/2,'DataScale',1));
        faileddown=downcrosses(~ismembertol(downcrosses,rewlocs,round(Fs(m,e))/2,'DataScale',1));
        leave_bottom_after_reward=[];
        leave_bottom_wo_reward=[];
        for i=1:length(downleaving)
            after_rew=downleaving(i)-rewardsdown;
            after_fail=downleaving(i)-faileddown;
            after_rew=after_rew(after_rew>0);
            after_fail=after_fail(after_fail>0);
            if isempty(after_fail)
                leave_bottom_after_reward=[leave_bottom_after_reward, downleaving(i)];
            elseif isempty(after_rew)
                leave_bottom_wo_reward=[leave_bottom_wo_reward,downleaving(i)];
            elseif min(min(after_fail),min(after_rew))==min(after_rew)
                leave_bottom_after_reward=[leave_bottom_after_reward, downleaving(i)];
            else leave_bottom_wo_reward=[leave_bottom_wo_reward,downleaving(i)];
            end
        end
        
        
        upcrosses=find((diff(y1(:,1)-170>0))==1);
        upleaving=find((diff(y1(:,1)-170>0))==-1);
        rewardsup=upcrosses(ismembertol(upcrosses,rewlocs,round(Fs(m,e))/2,'DataScale',1));
        failedup=upcrosses(~ismembertol(upcrosses,rewlocs,round(Fs(m,e))/2,'DataScale',1));
        leave_top_after_reward=[];
        leave_top_wo_reward=[];
        for i=1:length(upleaving)
            after_rew=upleaving(i)-rewardsup;
            after_fail=upleaving(i)-failedup;
            after_rew=after_rew(after_rew>0);
            after_fail=after_fail(after_fail>0);
            if isempty(after_fail)
                leave_top_after_reward=[leave_top_after_reward, upleaving(i)];
            elseif isempty(after_rew)
                leave_top_wo_reward=[leave_top_wo_reward,upleaving(i)];
            elseif min(min(after_fail),min(after_rew))==min(after_rew)
                leave_top_after_reward=[leave_top_after_reward, upleaving(i)];
            else leave_top_wo_reward=[leave_top_wo_reward,upleaving(i)];
            end
        end
        
        m
        e
        mrew_do_F=nan;mrew_do_V=nan;semrew_do_F=nan;semrew_do_V=nan;rew_do_vf=nan;rew_do_vr=nan;
        mfail_do_F=nan;mfail_do_V=nan;semfail_do_F=nan;semfail_do_V=nan;fail_do_vf=nan;fail_do_vr=nan;
        mrew_lb_F=nan;mrew_lb_V=nan;semrew_lb_F=nan;semrew_lb_V=nan;rew_lb_vf=nan;rew_lb_vr=nan;
        mworew_lb_F=nan;mworew_lb_V=nan;semworew_lb_F=nan;semworew_lb_V=nan;worew_lb_vf=nan;worew_lb_vr=nan;
        mrew_lt_F=nan;mrew_lt_V=nan;semrew_lt_F=nan;semrew_lt_V=nan;rew_lt_vf=nan;rew_lt_vr=nan;
        mworew_lt_F=nan;mworew_lt_V=nan;semworew_lt_F=nan;semworew_lt_V=nan;worew_lt_vf=nan;worew_lt_vr=nan;
        mrew_up_F=nan;mrew_up_V=nan;semrew_up_F=nan;semrew_up_V=nan;rew_up_vf=nan;rew_up_vr=nan;
        mfail_up_F=nan;mfail_up_V=nan;semfail_up_F=nan;semfail_up_V=nan;fail_up_vf=nan;fail_up_vr=nan;
        
        if ~isempty(rewardsdown)
            [rew_do_F,rew_do_V,mrew_do_F,mrew_do_V,semrew_do_F,semrew_do_V,rew_do_vf,rew_do_vr]=avg_activity(rewardsdown,Fs(m,e),F1,f1,r1);
        end
        if ~isempty(faileddown)
            [fail_do_F,fail_do_V,mfail_do_F,mfail_do_V,semfail_do_F,semfail_do_V,fail_do_vf,fail_do_vr]=avg_activity(faileddown,Fs(m,e),F1,f1,r1);
        end
        if ~isempty(leave_bottom_after_reward)
            [rew_lb_F,rew_lb_V,mrew_lb_F,mrew_lb_V,semrew_lb_F,semrew_lb_V,rew_lb_vf,rew_lb_vr]=avg_activity(leave_bottom_after_reward,Fs(m,e),F1,f1,r1);
        end
        if ~isempty(leave_bottom_wo_reward)
            [worew_lb_F,worew_lb_V,mworew_lb_F,mworew_lb_V,semworew_lb_F,semworew_lb_V,worew_lb_vf,worew_lb_vr]=avg_activity(leave_bottom_wo_reward,Fs(m,e),F1,f1,r1);
        end
        if ~isempty(rewardsup)
            [rew_up_F,rew_up_V,mrew_up_F,mrew_up_V,semrew_up_F,semrew_up_V,rew_up_vf,rew_up_vr]=avg_activity(rewardsup,Fs(m,e),F1,f1,r1);
        end
        if ~isempty(failedup)
            [fail_up_F,fail_up_V,mfail_up_F,mfail_up_V,semfail_up_F,semfail_up_V,fail_up_vf,fail_up_vr]=avg_activity(failedup,Fs(m,e),F1,f1,r1);
        end
        if ~isempty(leave_top_after_reward)
            [rew_lt_F,rew_lt_V,mrew_lt_F,mrew_lt_V,semrew_lt_F,semrew_lt_V,rew_lt_vf,rew_lt_vr]=avg_activity(leave_top_after_reward,Fs(m,e),F1,f1,r1);
        end
        if ~isempty(leave_top_wo_reward)
            [worew_lt_F,worew_lt_V,mworew_lt_F,mworew_lt_V,semworew_lt_F,semworew_lt_V,worew_lt_vf,worew_lt_vr]=avg_activity(leave_top_wo_reward,Fs(m,e),F1,f1,r1);
        end
        
        
        
        maxinfo=zeros(2,size(F1,2));
        timeleaving=linspace(-4,6+1/Fs(m,e),1+10*Fs(m,e));
        timeentering=linspace(-6,4+1/Fs(m,e),1+10*Fs(m,e));
        for cellnum=1:size(F1,2)
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            suptitle(['Cell: ',num2str(cellnum)])
            subplot(3,4,1)
            %             if ~isnan(mrew_do_F)
            plot(timeentering,mrew_do_F(:,cellnum))
            hold on
            plot(timeentering,mrew_do_F(:,cellnum)+1.96*(semrew_do_F(:,cellnum)),'g:')
            plot(timeentering,mrew_do_F(:,cellnum)-1.96*(semrew_do_F(:,cellnum)),'g:')
            line([0,0],[min(mrew_do_F(:,cellnum)-1.96*(semrew_do_F(:,cellnum))),max(mrew_do_F(:,cellnum)+1.96*(semrew_do_F(:,cellnum)))])
            title('dF Successfully entering bottom')
            %             end
            %             if ~isnan(mfail_do_F)
            subplot(3,4,2)
            if numel(mfail_do_F)>1 &&~isnan(nanmean(mfail_do_V(:,cellnum))) && ~isnan(nanmean(mfail_do_F(:,cellnum)))
                plot(timeentering,mfail_do_F(:,cellnum))
                hold on
                plot(timeentering,mfail_do_F(:,cellnum)+1.96*(semfail_do_F(:,cellnum)),'g:')
                plot(timeentering,mfail_do_F(:,cellnum)-1.96*(semfail_do_F(:,cellnum)),'g:')
                line([0,0],[min(mfail_do_F(:,cellnum)-1.96*(semfail_do_F(:,cellnum))),max(mrew_do_F(:,cellnum)+1.96*(semfail_do_F(:,cellnum)))])
            end
            title('dF Unsuccessfully entering bottom')
            %             end
            %             if ~isnan(mrew_do_F)
            subplot(3,4,5)
            plot(timeentering,mrew_do_V(:,cellnum))
            hold on
            plot(timeentering,mrew_do_V(:,cellnum)+1.96*(semrew_do_V(:,cellnum)),'g:')
            plot(timeentering,mrew_do_V(:,cellnum)-1.96*(semrew_do_V(:,cellnum)),'g:')
            plot(timeentering,rew_do_vf(:,cellnum))
            plot(timeentering,rew_do_vr(:,cellnum))
            line([0,0],[min(mrew_do_V(:,cellnum)-1.96*(semrew_do_V(:,cellnum))),max(mrew_do_V(:,cellnum)+1.96*(semrew_do_V(:,cellnum)))])
            title('Speed Successfully entering bottom')
            %             end
            %             if ~isnan(mfail_do_F)
            subplot(3,4,6)
            if ~isnan(nanmean(mfail_do_V))
                plot(timeentering,mfail_do_V(:,cellnum))
                hold on
                plot(timeentering,mfail_do_V(:,cellnum)+1.96*(semfail_do_V(:,cellnum)),'g:')
                plot(timeentering,mfail_do_V(:,cellnum)-1.96*(semfail_do_V(:,cellnum)),'g:')
                plot(timeentering,fail_do_vf(:,cellnum))
                plot(timeentering,fail_do_vr(:,cellnum))
                line([0,0],[min(mfail_do_V(:,cellnum)-1.96*(semfail_do_V(:,cellnum))),max(mfail_do_V(:,cellnum)+1.96*(semfail_do_V(:,cellnum)))])
            end
            title('Speed Unsuccessfully entering bottom')
            
            subplot(3,4,9)
            if numel(mrew_do_F)>1 && ~isnan(nanmean(mrew_do_V(:,cellnum))) && ~isnan(nanmean(mrew_do_F(:,cellnum)))
            [xc,lag]=xcorr(mrew_do_F(:,cellnum)-mean(mrew_do_F(:,cellnum)),mrew_do_V(:,cellnum)-mean(mrew_do_V(:,cellnum)),'coeff');
            plot(lag,xc)
            hold on
            set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
            line([0 0], [min(xc),max(xc)])
            xlabel('Seconds')
            [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
            maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
            scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
            text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            end
            title('XCorr Successfully entering bottom')
            %             end
            %             if ~isnan(mfail_do_F)
            subplot(3,4,10)
            if ~isnan(nanmean(mfail_do_V))
                [xc,lag]=xcorr(mfail_do_F(:,cellnum)-mean(mfail_do_F(:,cellnum)),mfail_do_V(:,cellnum)-mean(mfail_do_F(:,cellnum)),'coeff');
                plot(lag,xc)
                hold on
                set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
                line([0 0], [min(xc),max(xc)])
                xlabel('Seconds')
                [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
                maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
                scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
                text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            end
            title('XCorr Unsuccessfully entering bottom')
            %             if ~isnan(mrew_up_F)
            
            subplot(3,4,3)
            if numel(mrew_up_F)>1 &&~isnan(nanmean(mrew_up_V(:,cellnum))) && ~isnan(nanmean(mrew_up_F(:,cellnum)))
                plot(timeentering,mrew_up_F(:,cellnum))
                hold on
                plot(timeentering,mrew_up_F(:,cellnum)+1.96*(semrew_up_F(:,cellnum)),'g:')
                plot(timeentering,mrew_up_F(:,cellnum)-1.96*(semrew_up_F(:,cellnum)),'g:')
                line([0,0],[min(mrew_up_F(:,cellnum)-1.96*(semrew_up_F(:,cellnum))),max(mrew_up_F(:,cellnum)+1.96*(semrew_up_F(:,cellnum)))])
            end
            title('dF Successfully entering top')
            %             end
            subplot(3,4,4)
            if numel(mfail_up_F)>1 && ~isnan(nanmean(mfail_up_V(:,cellnum))) && ~isnan(nanmean(mfail_up_F(:,cellnum)))
                plot(timeentering,mfail_up_F(:,cellnum))
                hold on
                plot(timeentering,mfail_up_F(:,cellnum)+1.96*(semfail_up_F(:,cellnum)),'g:')
                plot(timeentering,mfail_up_F(:,cellnum)-1.96*(semfail_up_F(:,cellnum)),'g:')
                line([0,0],[min(mfail_up_F(:,cellnum)-1.96*(semfail_up_F(:,cellnum))),max(mrew_up_F(:,cellnum)+1.96*(semfail_up_F(:,cellnum)))])
            end
            title('dF Unsuccessfully entering top')
            subplot(3,4,7)
            if ~isnan(nanmean(mrew_up_V))
                
                plot(timeentering,mrew_up_V(:,cellnum))
                hold on
                plot(timeentering,mrew_up_V(:,cellnum)+1.96*(semrew_up_V(:,cellnum)),'g:')
                plot(timeentering,mrew_up_V(:,cellnum)-1.96*(semrew_up_V(:,cellnum)),'g:')
                plot(timeentering,rew_up_vf(:,cellnum))
                plot(timeentering,rew_up_vr(:,cellnum))
                line([0,0],[min(mrew_up_V(:,cellnum)-1.96*(semrew_up_V(:,cellnum))),max(mrew_up_V(:,cellnum)+1.96*(semrew_up_V(:,cellnum)))])
                title('Speed Successfully entering top')
            end
            subplot(3,4,8)
            if numel(mfail_up_V)>1 && ~isnan(nanmean(mfail_up_V(:,cellnum)))
                plot(timeentering,mfail_up_V(:,cellnum))
                hold on
                plot(timeentering,mfail_up_V(:,cellnum)+1.96*(semfail_up_V(:,cellnum)),'g:')
                plot(timeentering,mfail_up_V(:,cellnum)-1.96*(semfail_up_V(:,cellnum)),'g:')
                plot(timeentering,fail_up_vf(:,cellnum))
                plot(timeentering,fail_up_vr(:,cellnum))
                line([0,0],[min(mfail_up_V(:,cellnum)-1.96*(semfail_up_V(:,cellnum))),max(mfail_up_V(:,cellnum)+1.96*(semfail_up_V(:,cellnum)))])
            end
            title('Speed Unsuccessfully entering top')
            
            subplot(3,4,11)
            if ~isnan(nanmean(mrew_up_F))
                [xc,lag]=xcorr(mrew_up_F(:,cellnum)-mean(mrew_up_F(:,cellnum)),mrew_up_V(:,cellnum)-mean(mrew_up_V(:,cellnum)),'coeff');
                plot(lag,xc)
                hold on
                set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
                line([0 0], [min(xc),max(xc)])
                xlabel('Seconds')
                [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
                maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
                scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
                text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            end
            title('XCorr Successfully entering top')
            
            
            subplot(3,4,12)
            if numel(mfail_up_F)>1 && ~isnan(nanmean(mfail_up_V(:,cellnum))) && ~isnan(nanmean(mfail_up_F(:,cellnum)))
                [xc,lag]=xcorr(mfail_up_F(:,cellnum)-mean(mfail_up_F(:,cellnum)),mfail_up_V(:,cellnum)-mean(mfail_up_V(:,cellnum)),'coeff');
                plot(lag,xc)
                hold on
                set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
                line([0 0], [min(xc),max(xc)])
                xlabel('Seconds')
                [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
                maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
                scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
                text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            end
            title('XCorr Unsuccessfully entering top')
            
            
            saveas(gca,[saveDir,'\', 'Activity entering endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.jpg']);
            savefig([saveDir,'\', 'Activity entering endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.fig']);
            
            
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            suptitle(['Cell: ',num2str(cellnum)])
            subplot(2,4,1)
            %             if ~isnan(mrew_do_F)
            imagesc(squeeze(rew_do_F(:,cellnum,:))')
            t=textscan(num2str(-6:4),'%s');
            set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            colormap('jet')
            hold on
            title('dF Successfully entering bottom')
            %             end
            %             if ~isnan(mfail_do_F)
            subplot(2,4,2)
            if numel(mfail_do_F)>1 && ~isnan(nanmean(mfail_do_F(:,cellnum)))
                imagesc(squeeze(fail_do_F(:,cellnum,:))')
                t=textscan(num2str(-6:4),'%s');
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            end
            title('dF Unsuccessfully entering bottom')
            %             end
            %             if ~isnan(mrew_do_F)
            subplot(2,4,5)
            imagesc(squeeze(rew_do_V(:,cellnum,:))')
            t=textscan(num2str(-6:4),'%s');
            set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            colormap('jet')
            title('Speed Successfully entering bottom')
            %             end
            %             if ~isnan(mfail_do_F)
            subplot(2,4,6)
            if ~isnan(nanmean(mfail_do_V))
                imagesc(squeeze(fail_do_V(:,cellnum,:))')
                t=textscan(num2str(-6:4),'%s');
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
                colormap('jet')
            end
            title('Speed Unsuccessfully entering bottom')
            
            %             if ~isnan(mrew_up_F)
            subplot(2,4,3)
            if ~isnan(nanmean(mrew_up_V))
                imagesc(squeeze(rew_up_F(:,cellnum,:))')
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            end
            title('dF Successfully entering top')
            %             end
            subplot(2,4,4)
            if numel(mfail_up_F)>1 && ~isnan(nanmean(mfail_up_F(:,cellnum)))
                imagesc(squeeze(fail_up_F(:,cellnum,:))')
                t=textscan(num2str(-6:4),'%s');
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
                colormap('jet')
            end
            title('dF Unsuccessfully entering top')
            subplot(2,4,7)
            if ~isnan(nanmean(mrew_up_V))
                
                imagesc(squeeze(rew_up_V(:,cellnum,:))')
                t=textscan(num2str(-6:4),'%s');
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
                colormap('jet')
                title('Speed Successfully entering top')
            end
            subplot(2,4,8)
            if numel(mfail_up_F)>1 && ~isnan(nanmean(mfail_up_F(:,cellnum)))
                imagesc(squeeze(fail_up_V(:,cellnum,:))')
                t=textscan(num2str(-6:4),'%s');
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
                colormap('jet')
            end
            title('Speed Unsuccessfully entering top')
            saveas(gca,[saveDir,'\', 'All traces entering endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.jpg']);
            savefig([saveDir,'\', 'All traces entering endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.fig']);
            
            
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            %             if ~isnan(mrew_lb_F)
            suptitle(['Cell: ',num2str(cellnum)])
            subplot(3,4,1)
            plot(timeleaving,mrew_lb_F(:,cellnum))
            hold on
            plot(timeleaving,mrew_lb_F(:,cellnum)+1.96*(semrew_lb_F(:,cellnum)),'g:')
            plot(timeleaving,mrew_lb_F(:,cellnum)-1.96*(semrew_lb_F(:,cellnum)),'g:')
            line([0,0],[min(mrew_lb_F(:,cellnum)-1.96*(semrew_lb_F(:,cellnum))),max(mrew_lb_F(:,cellnum)+1.96*(semrew_lb_F(:,cellnum)))])
            title('dF Leaving bottom after reward')
            %             end
            subplot(3,4,2)
            if numel(mworew_lb_F)>1 && ~isnan(nanmean(mworew_lb_F(:,cellnum)))
                plot(timeleaving,mworew_lb_F(:,cellnum))
                hold on
                plot(timeleaving,mworew_lb_F(:,cellnum)+1.96*(semworew_lb_F(:,cellnum)),'g:')
                plot(timeleaving,mworew_lb_F(:,cellnum)-1.96*(semworew_lb_F(:,cellnum)),'g:')
                line([0,0],[min(mworew_lb_F(:,cellnum)-1.96*(semworew_lb_F(:,cellnum))),max(mworew_lb_F(:,cellnum)+1.96*(semworew_lb_F(:,cellnum)))])
            end
            title('dF Leaving bottom after failure')
            %             if ~isnan(mrew_lb_V)
            subplot(3,4,5)
            plot(timeleaving,mrew_lb_V(:,cellnum))
            hold on
            plot(timeleaving,mrew_lb_V(:,cellnum)+1.96*(semrew_lb_V(:,cellnum)),'g:')
            plot(timeleaving,mrew_lb_V(:,cellnum)-1.96*(semrew_lb_V(:,cellnum)),'g:')
            plot(timeleaving,rew_lb_vf(:,cellnum))
            plot(timeleaving,rew_lb_vr(:,cellnum))
            line([0,0],[min(mrew_lb_V(:,cellnum)-1.96*(semrew_lb_V(:,cellnum))),max(mrew_lb_V(:,cellnum)+1.96*(semrew_lb_V(:,cellnum)))])
            title('Speed Leaving bottom after reward')
            %             end
            subplot(3,4,6)
            if ~isnan(nanmean(mworew_lb_V))
                plot(timeleaving,mworew_lb_V(:,cellnum))
                hold on
                plot(timeleaving,mworew_lb_V(:,cellnum)+1.96*(semworew_lb_V(:,cellnum)),'g:')
                plot(timeleaving,mworew_lb_V(:,cellnum)-1.96*(semworew_lb_V(:,cellnum)),'g:')
                plot(timeleaving,worew_lb_vf(:,cellnum))
                plot(timeleaving,worew_lb_vr(:,cellnum))
                line([0,0],[min(mworew_lb_V(:,cellnum)-1.96*(semworew_lb_V(:,cellnum))),max(mworew_lb_V(:,cellnum)+1.96*(semworew_lb_V(:,cellnum)))])
            end
            title('Speed leaving bottom after failure')
            
            subplot(3,4,9)
            [xc,lag]=xcorr(mrew_lb_F(:,cellnum)-mean(mrew_lb_F(:,cellnum)),mrew_lb_V(:,cellnum)-mean(mrew_lb_V(:,cellnum)),'coeff');
            plot(lag,xc)
            hold on
            set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
            line([0 0], [min(xc),max(xc)])
            xlabel('Seconds')
            [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
            maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
            scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
            text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            title('XCorr Successfully entering bottom')
            %             end
            %             if ~isnan(mfail_do_F)
            subplot(3,4,10)
            if ~isnan(nanmean(mworew_lb_V))
                [xc,lag]=xcorr(mworew_lb_F(:,cellnum)-mean(mworew_lb_F(:,cellnum)),mworew_lb_V(:,cellnum)-mean(mworew_lb_V(:,cellnum)),'coeff');
                plot(lag,xc)
                hold on
                set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
                line([0 0], [min(xc),max(xc)])
                xlabel('Seconds')
                [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
                maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
                scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
                text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            end
            title('XCorr Unsuccessfully entering bottom')
            
            subplot(3,4,3)
            if numel(mrew_lt_F)>1 &&~isnan(nanmean(mrew_lt_V(:,cellnum))) && ~isnan(nanmean(mrew_lt_F(:,cellnum)))
                
                plot(timeleaving,mrew_lt_F(:,cellnum))
                hold on
                plot(timeleaving,mrew_lt_F(:,cellnum)+1.96*(semrew_lt_F(:,cellnum)),'g:')
                plot(timeleaving,mrew_lt_F(:,cellnum)-1.96*(semrew_lt_F(:,cellnum)),'g:')
                line([0,0],[min(mrew_lt_F(:,cellnum)-1.96*(semrew_lt_F(:,cellnum))),max(mrew_lt_F(:,cellnum)+1.96*(semrew_lt_F(:,cellnum)))])
                title('dF leaving top after reward')
            end
            subplot(3,4,4)
            if numel(mworew_lt_F)>1 &&~isnan(nanmean(mrew_lt_V(:,cellnum))) && ~isnan(nanmean(mworew_lt_F(:,cellnum)))
                plot(timeleaving,mworew_lt_F(:,cellnum))
                hold on
                plot(timeleaving,mworew_lt_F(:,cellnum)+1.96*(semworew_lt_F(:,cellnum)),'g:')
                plot(timeleaving,mworew_lt_F(:,cellnum)-1.96*(semworew_lt_F(:,cellnum)),'g:')
                line([0,0],[min(mworew_lt_F(:,cellnum)-1.96*(semworew_lt_F(:,cellnum))),max(mworew_lt_F(:,cellnum)+1.96*(semworew_lt_F(:,cellnum)))])
            end
            title('dF leaving top after failure')
            
            subplot(3,4,7)
            if ~isnan(nanmean(mrew_lt_V))
                plot(timeleaving,mrew_lt_V(:,cellnum))
                hold on
                plot(timeleaving,mrew_lt_V(:,cellnum)+1.96*(semrew_lt_V(:,cellnum)),'g:')
                plot(timeleaving,mrew_lt_V(:,cellnum)-1.96*(semrew_lt_V(:,cellnum)),'g:')
                plot(timeleaving,rew_lt_vf(:,cellnum))
                plot(timeleaving,rew_lt_vr(:,cellnum))
                line([0,0],[min(mrew_lt_V(:,cellnum)-1.96*(semrew_lt_V(:,cellnum))),max(mrew_lt_V(:,cellnum)+1.96*(semrew_lt_V(:,cellnum)))])
                title('Speed leaving top after reward')
            end
            subplot(3,4,8)
            if ~isnan(nanmean(mworew_lt_V))
                plot(timeleaving,mworew_lt_V(:,cellnum))
                hold on
                plot(timeleaving,mworew_lt_V(:,cellnum)+1.96*(semworew_lt_V(:,cellnum)),'g:')
                plot(timeleaving,mworew_lt_V(:,cellnum)-1.96*(semworew_lt_V(:,cellnum)),'g:')
                plot(timeleaving,worew_lt_vf(:,cellnum))
                plot(timeleaving,worew_lt_vr(:,cellnum))
                line([0,0],[min(mworew_lt_V(:,cellnum)-1.96*(semworew_lt_V(:,cellnum))),max(mworew_lt_V(:,cellnum)+1.96*(semworew_lt_V(:,cellnum)))])
            end
            title('Speed leaving top after failure')
            
            subplot(3,4,11)
            if ~isnan(nanmean(mrew_lt_F))
                [xc,lag]=xcorr(mrew_lt_F(:,cellnum)-mean(mrew_lt_F(:,cellnum)),mrew_lt_V(:,cellnum)-mean(mrew_lt_V(:,cellnum)),'coeff');
                plot(lag,xc)
                hold on
                set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
                line([0 0], [min(xc),max(xc)])
                xlabel('Seconds')
                [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
                maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
                scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
                text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            end
            title('XCorr Successfully entering top')
            
            
            subplot(3,4,12)
            if numel(mworew_lt_F)>1 && ~isnan(nanmean(mworew_lt_V(:,cellnum))) && ~isnan(nanmean(mworew_lt_F(:,cellnum)))
                [xc,lag]=xcorr(mworew_lt_F(:,cellnum)-mean(mworew_lt_F(:,cellnum)),mworew_lt_V(:,cellnum)-mean(mworew_lt_V(:,cellnum)),'coeff');
                plot(lag,xc)
                hold on
                set(gca,'XTick',lag(1):5*Fs(m,e):lag(end),'XTickLabel',strsplit(num2str((lag(1):5*Fs(m,e):lag(end))/(Fs(m,e)))))
                line([0 0], [min(xc),max(xc)])
                xlabel('Seconds')
                [maxinfo(1,cellnum),maxinfo(2,cellnum)]=max(xc);
                maxinfo(2,cellnum)=lag(maxinfo(2,cellnum));
                scatter(maxinfo(2,cellnum),maxinfo(1,cellnum),'k+')
                text(maxinfo(2,cellnum),maxinfo(1,cellnum)+.02*maxinfo(1,cellnum),num2str(maxinfo(2,cellnum)/Fs(m,e)))
            end
            title('XCorr Unsuccessfully entering top')
            
            saveas(gca,[saveDir,'\', 'Activity leaving endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.jpg']);
            savefig([saveDir,'\', 'Activity leaving endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.fig']);
            t=textscan(num2str(-4:6),'%s');
            
            
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            %             if ~isnan(mrew_lb_F)
            suptitle(['Cell: ',num2str(cellnum)])
            subplot(2,4,1)
            imagesc(squeeze(rew_lb_F(:,cellnum,:))')
            colormap('jet')
            set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            title('dF Leaving bottom after reward')
            %             end
            subplot(2,4,2)
            if numel(mworew_lb_F)>1 && ~isnan(nanmean(mworew_lb_V(:,cellnum))) && ~isnan(nanmean(mworew_lb_F(:,cellnum)))
                imagesc(squeeze(worew_lb_F(:,cellnum,:))')
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            end
            title('dF Leaving bottom after failure')
            %             if ~isnan(mrew_lb_V)
            subplot(2,4,5)
            imagesc(squeeze(rew_lb_V(:,cellnum,:))')
            set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            title('Speed Leaving bottom after reward')
            %             end
            subplot(2,4,6)
            if ~isnan(nanmean(mworew_lb_V))
                imagesc(squeeze(worew_lb_V(:,cellnum,:))')
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            end
            title('Speed leaving bottom after failure')
            subplot(2,4,3)
            if numel(mrew_lt_F)>1 &&~isnan(nanmean(mrew_lt_V(:,cellnum))) && ~isnan(nanmean(mrew_lt_F(:,cellnum)))
                imagesc(squeeze(rew_lt_F(:,cellnum,:))')
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
                title('dF leaving top after reward')
            end
            subplot(2,4,4)
            if numel(mworew_lt_F)>1 &&~isnan(nanmean(mworew_lt_V(:,cellnum))) && ~isnan(nanmean(mworew_lt_F(:,cellnum)))
                imagesc(squeeze(worew_lt_F(:,cellnum,:))')
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            end
            title('dF leaving top after failure')
            
            subplot(2,4,7)
            if numel(mrew_lt_V)>1 && ~isnan(nanmean(mrew_lt_V(:,cellnum)))
                imagesc(squeeze(rew_lt_V(:,cellnum,:))')
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
                title('Speed leaving top after reward')
            end
            subplot(2,4,8)
            if ~isnan(nanmean(mworew_lt_V))
                imagesc(squeeze(worew_lt_V(:,cellnum,:))')
                set(gca,'XTick',1:Fs(m,e):Fs(m,e)*10,'XTickLabels',t{1})
            end
            title('Speed leaving top after failure')
            saveas(gca,[saveDir,'\', 'All Traces leaving endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.jpg']);
            savefig([saveDir,'\', 'All Traces leaving endzone',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.fig']);
            
        end
        close all
    end
end