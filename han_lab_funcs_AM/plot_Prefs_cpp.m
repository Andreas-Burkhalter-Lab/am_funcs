function plot_Prefs_cpp(yb,up,down,Sp,pname,name,type,varargin)
if nargin>7
   part=varargin{1}; 
end

useinds=1:round(5.2*60*15);

%Bottom is A, Top is B
percentInB=sum(up)/(sum(up)+sum(down));
percentInA=sum(down)/(sum(up)+sum(down));
uA=mean(Sp(down,:))';
uB=mean(Sp(up,:))';
pref_A=(uA-uB)./(uA+uB);

startPup=find(diff(bwlabel(up))>0);
stopPup=find(diff(bwlabel(up))<0);
if length(startPup)<length(stopPup) || up(1)==1
startPup=[1; startPup];
end
if length(startPup)>length(stopPup)
stopPup=[stopPup; length(yb)];
end
startPdo=find(diff(bwlabel(down))>0);
stopPdo=find(diff(bwlabel(down))<0);

if length(startPdo)<length(stopPdo) || down(1)==1
startPdo=[1; startPdo];
end
if length(startPdo)>length(stopPdo)
stopPdo=[stopPdo; length(yb)];
end

%% Bar Graph With Residency Times
figure('units','normalized', 'Position', [.01 .05 .95 .5]);
subplot(1,3,1);
bar(1,percentInA,'r'); hold on
bar(2,percentInB,'b')
set(gca,'XTick',1:2,'XTickLabel',{'Env A','Env B'})
title('Mouse Residency Times')
ylabel('Percent Time in Region')
xlim([.5 2.5])
ylim([0 .9])
% saveas(gca,[pname,name,'\',name,' ',type,' Mouse Residency.jpg']);

%% Color-coded stem plot. Positive is A is red
[stA,indSpecSp]=sort(pref_A);
stB=stA;
stA(stA<0)=nan;
stB(stB>0)=nan;
% figure; stem(stA,'color','red','MarkerSize',.1);
% hold on; stem(stB,'color','b','MarkerSize',.1);
% xlabel('Cells')
% ylabel('Region Preference')

%% Bar Graph With Cell Percentages Times
subplot(1,3,2);
bar(1,sum(pref_A>0)/length(pref_A),'r'); hold on
bar(2,sum(pref_A<0)/length(pref_A),'b')
xlim([.5 2.5])
ylim([0 .8])

set(gca,'XTick',[1 2],'XTickLabel',{'Env A','Env B'})
title('Cell Preference Ratios')
ylabel('Percent Of Population Prefering Region')

%% Histogram with Cell Specificity
subplot(1,3,3);
h1=histogram(-pref_A(pref_A>0),-1:.1:1,'FaceColor','r');
hold on
h2=histogram(-pref_A(pref_A<=0),-1:.1:1,'FaceColor','b');
ylim([0 80])
title('Histogram of Cell Specificity Values')
xlabel('Cell Specificity Value')
ylabel('Number of Cells')
saveas(gca,[pname,name,'\',name,' ',type,' Cell Specificity Bars.jpg']);
saveas(gca,[pname,name,'\',name,' ',type,' Cell Specificity Bars.svg']);

%% 
figure('units','normalized', 'Position', [.01 .05 .95 .86]);
s1=subplot(2,10,1:9);
imagesc(Sp(useinds,indSpecSp)');
title('Cell Raster Organized by Region Specificity')
colormap(flipud(gray))
for ii=1:length(startPup)
patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 size(Sp,2) size(Sp,2) 0],...
'b','FaceAlpha',.1,'EdgeAlpha',0)
end
for ii=1:length(startPdo)
patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 size(Sp,2) size(Sp,2) 0],...
'r','FaceAlpha',.1,'EdgeAlpha',0)
end
ylabel('Cells')
subplot(2,10,10);
stem(stA,'color','red','MarkerSize',.1);
title('Region Specificity')
hold on; stem(stB,'color','b','MarkerSize',.1);
view(90,90)
xlim([1 length(stA)])
s2=subplot(2,10,11:19);
plot(yb,'k')
title('Mouse Position')
ylabel('Region')
axis(gca,'tight')
xlabel('Frames')
set(gca,'YTick',[1 2],'YTickLabel',{'A','B'})
for ii=1:length(startPup)
patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 3 3 0],...
'b','FaceAlpha',.1,'EdgeAlpha',0)
end
for ii=1:length(startPdo)
patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 3 3 0],...
'r','FaceAlpha',.1,'EdgeAlpha',0)
end
linkaxes([s1 s2],'x')

saveas(gca,[pname,name,'\',name,' ',type,' RasterWithRegions.jpg']);
saveas(gca,[pname,name,'\',name,' ',type,' RasterWithRegions.svg']);


if nargin>7
    import mlreportgen.dom.*;
    I=Image([pname,name,'\',name,' ',type,' Cell Specificity Bars.jpg']);
    append(part,TableEntry(I));  
    I=Image([pname,name,'\',name,' ',type,' RasterWithRegions.jpg']);
    append(part,TableEntry(I));
end



end