function[]=plotcompilemsac(savefilename,xlimmaxes);

% function[]=plotcompilemac(savefilename,xlimmaxes);
% Function to plot the previously saved results of compilemsac.m, which
%	compiles microsaccade statistics from a long list of experiments
%	for an animal.  Savefilename is the name of the file that
%	compilemsac.m saves; xlimmaxes is a three element optional input
%	that tells program xlimit maxes to use instead of defaults.
% LOADS:  file in /HOME/crista/files/datfiles/sac
% PLOTS:  Plots distributions of saccade durations, saccade vector 
%	amplitudes, and intersaccade intervals (real and pseudo).

set(0,'ShowHiddenHandles','off')

cd /HOME/crista/files/datfiles/sac
eval(['load ',savefilename])

if nargin<2;xlimmaxes=[];end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plotting Figure %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','inches','Position',[9 .3 5 9])
set(gcf,'DefaultAxesFontName','Palatino','DefaultAxesFontSize',8)
set(gcf,'DefaultTextFontName','Palatino','DefaultTextFontSize',10)
set(gcf,'DefaultAxesNextPlot','add','DefaultAxesUnits','inches')
set(gcf,'DefaultAxesTickDir','out')
set(gcf,'PaperUnits','inches','PaperPosition',[1.75 2.25 8 10.5])

axes('Position',[.5 7.8 4 .5])
axes('Position',[.5 5.5 4 2])
axes('Position',[.5 3 4 2])
axes('Position',[.5 .5 4 2])
c=get(gcf,'children');c=sort(c);
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Getting Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Getting Saccade Rate
msaccount=size(msacstatsall,1);
tottime=sum(trlstatsall(:,3)-trlstatsall(:,2));
msacrate=msaccount/(tottime/1000);

sacdurs=msacstatsall(:,3);
sacamps=msacstatsall(:,4);
rimsis=intstatsall(find(intstatsall(:,3)==1),1);
pimsis=intstatsall(find(intstatsall(:,3)==0),1);
aimsis=[rimsis;pimsis];

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plotting Data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
	
axes(c(1))
text('Position',[1 1],'String',['saccade count:  ',num2str(msaccount)],'HorizontalAlignment','right')
text('Position',[1 .6],'String',['saccade rate:  ',num2str(msacrate),' sac/sec'],'HorizontalAlignment','right')
text('Position',[0 1],'String',['Saccade Statistics'],'FontSize',12)
text('Position',[0 .6],'String',['File:  ',savefilename],'FontSize',10)
text('Position',[0 .2],'String',['prog:  plotcompilemsac.m'],'FontSize',8)
axis('off')
	
axes(c(2))
maxbin=ceil(max(sacdurs)/10)*10;
if ~isempty(xlimmaxes);maxbin=xlimmaxes(1);end
hist(sacdurs,[maxbin/20:maxbin/20:maxbin])
set(gca,'XLim',[0 maxbin*1.05],'XTick',[0:maxbin/10:maxbin])
d=get(gca,'children');set(d(1),'FaceColor',[.7 .7 .7]);
title('saccade durations (msec)','FontSize',10)

axes(c(3))
maxbin=ceil(max(sacamps));
if ~isempty(xlimmaxes);maxbin=xlimmaxes(2);end
hist(sacamps,[maxbin/20:maxbin/20:maxbin])
set(gca,'XLim',[0 maxbin*1.05],'XTick',[0:maxbin/10:maxbin])
d=get(gca,'children');set(d(1),'FaceColor',[.7 .7 .7]);
title('saccade vector amplitudes (degrees)','FontSize',10)
	
axes(c(4))
maxbin=ceil(max(aimsis)/100)*100;
if ~isempty(xlimmaxes);maxbin=xlimmaxes(3);end
abins=[maxbin/20:maxbin/20:maxbin];
[n]=hist(rimsis,abins);
[n2]=hist(pimsis,abins);
[n3]=hist(aimsis,abins);
bar(abins,[n;n2;n3]')
set(gca,'XLim',[0 maxbin*1.05],'XTick',[0:maxbin/10:maxbin])
colormap(gray)
title('intersaccade intervals (real and pseudo; msec)','FontSize',10)

