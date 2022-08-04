function[]=plotmsacalg(emh,emv,emt,trl,algparams)

% function[]=plotmsacalg(emh,emv,emt,trl,algparams)

%%%%% Bookkeeping %%%%%
cd /HOME/crista/matlab/eyes

%%%%% Checking Inputs %%%%%
if nargin<2;trl=1;end
if isempty(trl);trl=1;end
if nargin<3;algparams=[];end

%%%%% Getting Data %%%%%
[msac,emxout,demx,demh,demv,cemx,cemh,cemv,cemt,femx,femh,femv]=msacalg(emh,emv,emt,algparams);

algnum=algparams(1,1);
thresh=algparams(1,2);
tbmax=algparams(1,3);
tamin=algparams(1,4);
convfilt1=algparams(2,:);
if ~isempty(find(convfilt1))
	while convfilt1(length(convfilt1))==0;
		convfilt1=convfilt1(1:length(convfilt1)-1);
	end
else
	convfilt1=[];
end
if size(algparams,1)>=3
	convfilt2=algparams(3,:);
	while convfilt2(length(convfilt2))==0;
		convfilt2=convfilt2(1:length(convfilt2)-1);
	end
else
	convfilt2=[];
end

tin=find(~isnan(emt(trl,:)));
if ~isempty(cemh);cin=find(~isnan(cemh(trl,:)));
elseif ~isempty(cemx);cin=find(~isnan(cemx(trl,:)));end
if ~isempty(femh);fin=find(~isnan(femh(trl,:)));end
if ~isempty(demx);din=find(~isnan(demx(trl,:)));end
zin=find(~isnan(emxout(trl,:)));

%%%%% Making Figure %%%%%
figure('Units','inches','Position',[4.5 .1 8 10.5])
set(gcf,'DefaultAxesFontName','Palatino','DefaultAxesFontSize',7)
set(gcf,'DefaultTextFontName','Palatino','DefaultTextFontSize',8)
set(gcf,'DefaultAxesNextPlot','add','DefaultAxesUnits','inches')
set(gcf,'DefaultAxesTickDir','out','DefaultAxesTickLength',[.005 .025])
set(gcf,'PaperUnits','inches','PaperPosition',[.25 .25 8 10.5])
set(0,'ShowHiddenHandles','on')

axes('Position',[.5 5.9 7.1 3.3])
axes('Position',[.5 3.2 7.1 2.2])
axes('Position',[.5 0.5 7.1 2.2])
axes('Position',[.5 9.7 7.1 0.5])
c=get(gcf,'children');c=sort(c);

%%%%% Plotting Data %%%%%
msactimes=emt(trl,find(msac(trl,:)==1));
msacnum=length(msactimes);

axes(c(1))
plot(emt(trl,tin),emv(trl,tin),'m')
plot(emt(trl,tin),emh(trl,tin),'b')
ylabel('degrees')
set(gca,'XLim',[min(emt(trl,tin)) max(emt(trl,tin))])
set(gca,'XColor','w')
datmin=min([min(emv(trl,tin)) min(emh(trl,tin))]);
datmax=max([max(emv(trl,tin)) max(emh(trl,tin))]);
offset=(datmax+datmin)/2;
offset=round(offset*10)/10;
dathalf=(datmax-datmin)/2;
dathalf=ceil(dathalf*10)/10;
if dathalf<.5;dathalf=.5;end
ymin=offset-dathalf;
ymax=offset+dathalf;
set(gca,'YLim',[ymin ymax],'YTick',[ymin:.2:ymax])	
title('emh,emv')

if (algnum==2 | algnum==3)
	plot(emt(trl,fin),femv(trl,fin)+.4,'m')
	plot(emt(trl,fin),femh(trl,fin)+.4,'b')
	set(gca,'YLim',[ymin ymax+.5],'YTick',[ymin:.2:ymax+.6])
	title('emh,emv and femh,femv')	
end

axes(c(2))
if (algnum==1 | algnum==3)
	plot(emt(trl,cin),cemv(trl,cin),'m')
	plot(emt(trl,cin),cemh(trl,cin),'b')
	title('cemh,cemv')
elseif algnum==2
	plot(emt(trl,cin),cemx(trl,cin),'k')
	title('cemx (diff of sqrt(femh.^2+femv.^2))')
end
set(gca,'XLim',[min(emt(trl,tin)) max(emt(trl,tin))])
set(gca,'XColor','w')

axes(c(3))
plot(emt(trl,zin),emxout(trl,zin),'k') 
set(gca,'XLim',[min(emt(trl,tin)) max(emt(trl,tin))])
set(gca,'YLim',[min(emxout(trl,zin)) max(emxout(trl,zin))])
xlabel('msec')
if algnum==1
	title('sqrt(cemh^2+cemv^2))')
elseif algnum==2
	title('abs(cemx/cemt)')
end

for ax=1:3
	axes(c(ax))
	ylimits=get(gca,'YLim');
	ymin=ylimits(1);ymax=ylimits(2);
	plot([msactimes;msactimes],[ymin*ones(1,msacnum);ymax*ones(1,msacnum)],'-','Color',[.7 .7 .7])
	d=get(gca,'children');d=sort(d);set(gca,'children',d);
end

axes(c(length(c)))
text('Position',[0 1],'String','Microsaccade Algorithm','FontSize',12)
text('Position',[0 .7],'String','program:  plotmsacalg.m','FontSize',7)
%text('Position',[.4 .9],'String',['file:    ',filename],'FontSize',10,'HorizontalAlignment','right')
text('Position',[.4 .9],'String',['algnum:    ',num2str(algnum)],'FontSize',10,'HorizontalAlignment','right')
text('Position',[.6 .9],'String',['thresh:    ',num2str(thresh)],'FontSize',10,'HorizontalAlignment','right')
text('Position',[.8 .9],'String',['tbmax:    ',num2str(tbmax)],'FontSize',10,'HorizontalAlignment','right')
text('Position',[1 .9],'String',['tamin:    ',num2str(tamin)],'FontSize',10,'HorizontalAlignment','right')
text('Position',[1 .4],'String',['convfilt1:    [',num2str(convfilt1),']'],'FontSize',10,'HorizontalAlignment','right')
text('Position',[1 0],'String',['convfilt2:    [',num2str(convfilt2),']'],'FontSize',10,'HorizontalAlignment','right')
axis('off')