function CPP_Pop(mouseCPP,ms,currs,pname,name,~)
anglePop{4}=[];
perInAPop{4}=[];
perInBPop{4}=[];
perPrefAPop{4}=[];
perPrefBPop{4}=[];
AllPrefA{4}=[];
AllPrefB{4}=[];
for m=ms
    for day=1:4
        
        temp=currs{m};
        e=temp(day);
        if ~isnan(e)
        %%
        maxh=876.2818;
        useinds=1:round(5.2*60*15);
        yb=(((mouseCPP(m).ybinned{e}{1}(useinds,1))/maxh)+1)*1.5;
        down=yb<1;
        up=yb>2;
        Sp=mouseCPP(m).spks{e}{1}(useinds,:);
        % Fc2=mouseCPP(m).FPyrs{curr}{1};
        F=calc_Fc3(mouseCPP(m).FPyrs{e}{1}(useinds,:));
        forward=mouseCPP(m).Forwards{e}{1}(useinds,:);
        rotation=mouseCPP(m).Rotations{e}{1}(useinds,:);
        %%
        %Bottom is A, Top is B
        percentInB=sum(up)/(sum(up)+sum(down));
        percentInA=sum(down)/(sum(up)+sum(down));
        uA=mean(F(down,:))';
        uB=mean(F(up,:))';
        pref_A=(uA-uB)./(uA+uB);
        %%
        perInAPop{day}=[perInAPop{day} percentInA];
        perInBPop{day}=[perInBPop{day} percentInB];
        %%
        perPrefAPop{day}=[perPrefAPop{day} sum(pref_A>0)/length(pref_A)];
        perPrefBPop{day}=[perPrefBPop{day} sum(pref_A<0)/length(pref_A)];
        %%
        AllPrefA{day}=[AllPrefA{day}; -pref_A(pref_A>0)];
        AllPrefB{day}=[AllPrefB{day}; -pref_A(pref_A<0)];
        %%
        Fup=F(up,:);
        Fdo=F(down,:);
        mFup=nanmean(Fup);
        mFdo=nanmean(Fdo);
        nanloc=isnan(mFup)|isnan(mFdo);
        mFup(nanloc)=[];
        mFdo(nanloc)=[];
        angle=dot(mFdo,mFup)/(norm(mFdo)*norm(mFup));
        anglePop{day}=[anglePop{day} angle];
        disp(['Mouse ',num2str(m),' Curr ', num2str(e), ' Time in A ',num2str(percentInA),...
            ' Time in B ',num2str(percentInB),' Percent Pref A ', num2str(sum(pref_A>0)/length(pref_A)),...
            ' Percent Pref B ',num2str(sum(pref_A<0)/length(pref_A))])
        end
    end
end
dayNums={'1','5','10','Test'};
figure('units','normalized', 'Position', [.01 .05 .9 .86]);
for day=1:4
    subplot(4,4,day)
    bar(1,mean(perInAPop{day}),'r'); hold on;
    bar(2,mean(perInBPop{day}),'b');
    A=mean_and_sem(perInAPop{day});
    B=mean_and_sem(perInBPop{day});
    errorbar([1 2],[A(1) B(1)],[A(2) B(2)],'.k');
    plot([perInAPop{day};perInBPop{day}],'k')
    xlim([.5 2.5])
    set(gca,'XTick',[1 2],'XTickLabel',{'Env A','Env B'})
%     title([name ,' ',dayNums{day}, '/n Mouse Residency Times'])
    title(sprintf([name ,' ',dayNums{day}, '\n Mouse Residency Times']))

    ylabel('Time in Region')
    
    subplot(4,4,day+4)
    bar(1,mean(perPrefAPop{day}),'r'); hold on;
    bar(2,mean(perPrefBPop{day}),'b');
    A=mean_and_sem(perPrefAPop{day});
    B=mean_and_sem(perPrefBPop{day});
        plot([perPrefAPop{day};perPrefBPop{day}],'k')

    errorbar([1 2],[A(1) B(1)],[A(2) B(2)],'.k');
    xlim([.5 2.5])
    ylim([0 1])
    set(gca,'XTick',[1 2],'XTickLabel',{'Env A','Env B'})
    title('Cell Preference Ratios')
    ylabel('Population Preference')
    
    subplot(4,4,day+8)
    hold on;
    h1=histogram(AllPrefA{day},-1:.1:1,'FaceColor','r');
    h2=histogram(AllPrefB{day},-1:.1:1,'FaceColor','b');
    ylim([0 200])
    title('Histogram of Cell Specificity Values')
    xlabel('Cell Specificity Value')
    ylabel('Number of Cells')
    
    subplot(4,4,day+12)
    line([0 1],[0 0],'color','k','linewidth',2);
    hold on;
    for ii=1:length(anglePop{day})
        line([0 cos(anglePop{day}(ii))], [0 sin(anglePop{day}(ii))],'color','k')
    end
    line([0 cos(mean(anglePop{day}))], [0 sin(mean(anglePop{day}))],'color','k','linewidth',2)
    ylim([-.1 1.1]); xlim([-.1 1.1])
    title(['Mean Angle: ',num2str(mean(anglePop{day})*180/pi)])
    xlabel([name ,' ',dayNums{day}])
    saveas(gca,[pname,name,' Population Summary.jpg']);
        saveas(gca,[pname,name,' Population Summary.svg']);

end