function plot_CPP_across_days(mouseCPP,m,curr1,curr2,name)
pname='G:\MA Data\CPPs\Figures\';
maxh=876.2818;
useinds=1:round(5.2*60*15);
if ~exist([pname,'\',name,'\'],'dir')
    mkdir([pname,'\',name,'\']);
end

yb1=(((mouseCPP(m).ybinned{curr1}{1}(useinds,1))/maxh)+1)*1.5;
up1=yb1>2; down1=yb1<1;
yb2=(((mouseCPP(m).ybinned{curr2}{1}(useinds,1))/maxh)+1)*1.5;
up2=yb2>2; down2=yb2<1;
F1=calc_Fc3(mouseCPP(m).FPyrs{curr1}{1});
F2=calc_Fc3(mouseCPP(m).FPyrs{curr2}{1});

uA1=mean(F1(down1,:))';
uB1=mean(F1(up1,:))';
prefA1=(uA1-uB1)./(uA1+uB1);
uA2=mean(F1(down2,:))';
uB2=mean(F1(up2,:))';
prefA2=(uA2-uB2)./(uA2+uB2);
diffPref=prefA2-prefA1;
%%
figure('units','normalized', 'Position', [.01 .05 .95 .86]);
subplot(2,2,1)
bar(1,sum(down1)/(sum(up1)+sum(down1)),'r'); hold on
bar(2,sum(up1)/(sum(up1)+sum(down1)),'b')
set(gca,'XTick',1:2,'XTickLabel',{'Env A','Env B'})
title('Mouse Residency Times Day 1')
ylabel('Percent Time in Region')
xlim([.5 2.5])

subplot(2,2,2)
bar(1,sum(down2)/(sum(up2)+sum(down2)),'r'); hold on
bar(2,sum(up2)/(sum(up2)+sum(down2)),'b')
set(gca,'XTick',1:2,'XTickLabel',{'Env A','Env B'})
title('Mouse Residency Times Day 2')
ylabel('Percent Time in Region')
xlim([.5 2.5])

subplot(2,2,3)
bar(1,sum(prefA1>0)/length(prefA1),'r'); hold on
bar(2,sum(prefA1<0)/length(prefA1),'b')
set(gca,'XTick',1:2,'XTickLabel',{'Env A','Env B'})
title('Mouse Specificity Ratio Day 1')
ylabel('Percent Time in Region')
xlim([.5 2.5])

subplot(2,2,4)
bar(1,sum(prefA2>0)/length(prefA2),'r'); hold on
bar(2,sum(prefA2<0)/length(prefA2),'b')
set(gca,'XTick',1:2,'XTickLabel',{'Env A','Env B'})
title('Mouse Specificity Ratio Day 2')
ylabel('Percent Time in Region')
xlim([.5 2.5])




saveas(gca,[pname,name,'\',name,' Compare Mouse Ratios.jpg']);
%%
figure('units','normalized', 'Position', [.01 .05 .95 .86])
subplot(2,2,1)
[st1,indSpecSp1]=sort(prefA1);
stA1=st1; stB1=stA1;
stA1(stA1<0)=nan;
stB1(stB1>0)=nan;
stem(stA1,'color','red','MarkerSize',.1);
hold on; stem(stB1,'color','b','MarkerSize',.1);
xlim([1 length(st1)]); ylim([-1 1]);
xlabel('Cells')
ylabel('Region Preference')
title('Cell Preference on Day 1')
legend('Env A','Env B','Location','NW')

subplot(2,2,3)
[st2,~]=sort(prefA2);
stA2=st2; stB2=stA2;
stA2(stA2<0)=nan;
stB2(stB2>0)=nan;
stem(stA2,'color','red','MarkerSize',.1);
hold on; stem(stB2,'color','b','MarkerSize',.1);
xlim([1 length(st1)]); ylim([-1 1]);
xlabel('Cells')
ylabel('Region Preference')
title('Cell Preference on Day 2')
legend('Env A','Env B','Location','NW')

subplot(2,2,2)
scatter(prefA1,prefA2,'k')
r=corrcoef(prefA1,prefA2);
xlabel('Day 1 Pref')
ylabel('Day 2 Pref')
title(['Cell Specificity Across Days: ',num2str(r(2,1))])


subplot(2,2,4)
reord=prefA2(indSpecSp1);
reordA=reord; reordB=reord;
reordA(isnan(stA1))=nan;
reordB(isnan(stB1))=nan;
stem(reordA,'color','red','MarkerSize',.1);
hold on; stem(reordB,'color','blue','MarkerSize',.1);
xlabel('Cells')
xlim([1 length(st1)]); ylim([-1 1]);
ylabel('Region Preference')
title('Day 2 Preference Ordered by Day 1')
legend('Env A','Env B','Location','NW')

saveas(gca,[pname,name,'\',name,' Compare Mouse Stems.jpg']);

%%
startPup=find(diff(bwlabel(up1))>0);
stopPup=find(diff(bwlabel(up1))<0);
if length(startPup)<length(stopPup)
startPup=[1; startPup];
elseif length(startPup)>length(stopPup)
stopPup=[stopPup; length(yb1)];
end
startPdo=find(diff(bwlabel(down1))>0);
stopPdo=find(diff(bwlabel(down1))<0);
if length(startPdo)<length(stopPdo)
startPdo=[1; startPdo];
elseif length(startPdo)>length(stopPdo)
stopPdo=[stopPdo; length(yb1)];
end

figure('units','normalized', 'Position', [.01 .05 .95 .86]);
s1=subplot(2,10,1:9);
imagesc(F1(useinds,indSpecSp1)')
colormap(flipud(gray))
for ii=1:length(startPup)
patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 size(F1,2) size(F1,2) 0],...
'b','FaceAlpha',.1,'EdgeAlpha',0)
end
for ii=1:length(startPdo)
patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 size(F1,2) size(F1,2) 0],...
'r','FaceAlpha',.1,'EdgeAlpha',0)
end
subplot(2,10,10);
stem(stA1,'color','red','MarkerSize',.1);
hold on; stem(stB1,'color','b','MarkerSize',.1);
view(90,90)
xlim([1 length(stA1)])
s2=subplot(2,10,11:19)
plot(yb1,'k')
axis(gca,'tight')
for ii=1:length(startPup)
patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 3 3 0],...
'b','FaceAlpha',.1,'EdgeAlpha',0)
end
for ii=1:length(startPdo)
patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 3 3 0],...
'r','FaceAlpha',.1,'EdgeAlpha',0)
end
linkaxes([s1 s2],'x')
saveas(gca,[pname,name,'\',name,'RasterWithRegionsDay1.jpg']);

%%
startPup=find(diff(bwlabel(up2))>0);
stopPup=find(diff(bwlabel(up2))<0);
if length(startPup)<length(stopPup)
startPup=[1; startPup];
elseif length(startPup)>length(stopPup)
stopPup=[stopPup; length(yb2)];
end
startPdo=find(diff(bwlabel(down2))>0);
stopPdo=find(diff(bwlabel(down2))<0);
if length(startPdo)<length(stopPdo)
startPdo=[1; startPdo];
elseif length(startPdo)>length(stopPdo)
stopPdo=[stopPdo; length(yb2)];
end
figure('units','normalized', 'Position', [.01 .05 .95 .86]);
s1=subplot(2,10,1:9);
imagesc(F2(useinds,indSpecSp1)')
colormap(flipud(gray))
for ii=1:length(startPup)
patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 size(F1,2) size(F1,2) 0],...
'b','FaceAlpha',.1,'EdgeAlpha',0)
end
for ii=1:length(startPdo)
patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 size(F1,2) size(F1,2) 0],...
'r','FaceAlpha',.1,'EdgeAlpha',0)
end
subplot(2,10,10);
stem(reordA,'color','red','MarkerSize',.1);
hold on; stem(reordB,'color','blue','MarkerSize',.1);
view(90,90)
xlim([1 length(stA1)])
s2=subplot(2,10,11:19)
plot(yb2,'k')
axis(gca,'tight')
for ii=1:length(startPup)
patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 3 3 0],...
'b','FaceAlpha',.1,'EdgeAlpha',0)
end
for ii=1:length(startPdo)
patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 3 3 0],...
'r','FaceAlpha',.1,'EdgeAlpha',0)
end
linkaxes([s1 s2],'x')
saveas(gca,[pname,name,'\',name,'RasterWithRegionsDay2OrderedByDay1.jpg']);















end