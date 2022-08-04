function collect_means=corr_EZactivity(mouse,varargin)
if nargin>1
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>4
        saveDir=varargin{4};
    else     saveDir='F:\MA Data\Interneurons\PubFigs\EZactivity\';
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
for m=mouseNums
    for e=expnum{m}
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
                leave_bottom_after_reward=[leave_after_reward, downleaving(i)];
            elseif isempty(after_rew)
                leave_bottom_wo_reward=[leave_bottom_wo_reward,downleaving(i)];
            elseif min(min(after_fail),min(after_rew))==min(after_rew)
                leave_bottom_after_reward=[leave_after_reward, downleaving(i)];
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
        
        
        [rew_do_F,rew_do_V,mrew_do_F,mrew_do_V,semrew_do_F,semrew_do_V]=avg_activity(rewardsdown,Fs(m,e),F1)
        [fail_do_F,fail_do_V,mfail_do_F,mfail_do_V,semfail_do_F,semfail_do_V]=avg_activity(faileddown,Fs(m,e),F1)
        [rew_lb_F,rew_lb_V,mrew_lb_F,mrew_lb_V,semrew_lb_F,semrew_lb_V]=avg_activity(leave_bottom_after_reward,Fs(m,e),F1)
        [worew_lb_F,worew_lb_V,mworew_lb_F,mworew_lb_V,semworew_lb_F,semworew_lb_V]=avg_activity(leave_bottom_wo_reward,Fs(m,e),F1)
        
        [rew_up_F,rew_up_V,mrew_up_F,mrew_up_V,semrew_up_F,semrew_up_V]=avg_activity(rewardsup,Fs(m,e),F1)
        [fail_up_F,fail_up_V,mfail_up_F,mfail_up_V,semfail_up_F,semfail_up_V]=avg_activity(failedup,Fs(m,e),F1)
        [rew_lt_F,rew_lt_V,mrew_lt_F,mrew_lt_V,semrew_lt_F,semrew_lt_V]=avg_activity(leave_top_after_reward,Fs(m,e),F1)
        [worew_lt_F,worew_lt_V,mworew_lt_F,mworew_lt_V,semworew_lt_F,semworew_lt_V]=avg_activity(leave_top_wo_reward,Fs(m,e),F1)
        
        
        
        
        
        timeentering=linspace(-1,5+1/Fs(m,e),binsize+1+Fs(m,e));
        timeleaving=linspace(-5,1+1/Fs(m,e),binsize+1+Fs(m,e));
        for cellnum=1:size(F1,2)
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        subplot(4,4,1)
        plot(timeentering,mrew_do_F(:,cellnum))
        hold on
        plot(timeentering,mrew_do_F(:,cellnum)+1.96(semrew_do_F(:,cellnum))
        title('dF Successfully entering bottom')
        
        saveas(gca,[saveDir,'\', 'Activity around bottom',num2str(m),',',num2str(e),' Cell ',num2str(cellnum),env,'.jpg']);
        savefig([saveDir,'\', 'XCorr from stop',num2str(m),',',num2str(e), ' All Cells ',env,'.fig']);
        end
    end
end

    function [F,V,mF,mV,semF,semV]=avg_activity(indices,Fs)
        binsize=4*Fs;
        F=zeros(binsize+Fs+1,size(F1,2),length(indices));
        V=zeros(binsize+Fs+1,size(F1,2),length(indices));
        for i=1:length(indices)
            start=max(1,indices(i)-binsize);
            stop=min(indices(i)+Fs,length(F1));
            tempF=F1(start:stop,:);
            tempV=sqrt(f1(start:stop,:).^2+r1(start:stop,:).^2);
            if (indices(i)-binsize)<1
                tempF=[nan(length(F)-length(tempF),size(F1,2));tempF];
                tempV=[nan(length(V)-length(tempV),size(F1,2));tempV];
            elseif (indices(i)+Fs)>length(F1)
                tempF=[tempF;nan(length(F)-length(tempF),size(F1,2))];
                tempV=[tempV;nan(length(V)-length(tempV),size(F1,2))];
            end
            F(:,:,i)=tempF;
            V(:,:,i)=tempV;
        end
        mF=mean(F,3);
        semF=std(F,0,3)/sqrt(size(F,3));
        mV=mean(F,3);
        semV=std(V,0,3)/sqrt(size(V,3));
    end