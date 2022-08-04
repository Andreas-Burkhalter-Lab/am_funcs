function[]=plotcompilemsac2(savefilename,xlimmaxes,trldurrange);

% function[]=plotcompilemac2(savefilename,xlimmaxes,trldurrange);
% Function to plot the previously saved results of compilemsac.m, which
%	compiles microsaccade statistics from a long list of experiments
%	for an animal.  Savefilename is the name of the file that
%	compilemsac.m saves; xlimmaxes is a three element optional input
%	that tells program xlimit maxes to use instead of defaults.
% LOADS:  file in /HOME/crista/files/datfiles/sac
% PLOTS:  Plots distributions of saccade durations, saccade vector 
%	amplitudes, and intersaccade intervals (real and pseudo).

cd /HOME/crista/files/datfiles/sac
eval(['load ',savefilename])

if nargin<2;xlimmaxes=[];end
if nargin<3;trldurrange=[];end

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
sacsts=msacstatsall(:,1);

if ~isempty(trldurrange)
	
	%%% last column of these variables is the file number index
	
	%%% column 1 is unique file/trl number, col 2 is trial dur in msec
	lastcol=size(trlstatsall,2);
	trldurs=[trlstatsall(:,lastcol)*1000+trlstatsall(:,1) trlstatsall(:,3)-trlstatsall(:,2)];
	
	%%% column 1 is unique file/trl number, col 2-5 is col 1-4 of msacstatsall
	lastcol=size(msacstatsall,2);
	msacstatsall2=[msacstatsall(:,lastcol)*1000+msacstatsall(:,5) msacstatsall(:,1:4)];
	
	%%% column 1 is unique file/trl number, col 2 is imsi, col 3 is rpflag
	lastcol=size(intstatsall,2);
	intstatsall2=[intstatsall(:,lastcol)*1000+intstatsall(:,2) intstatsall(:,[1 3])];
	
	trlstokeep=find(trldurs(:,2)>=trldurrange(1) & trldurs(:,2)<=trldurrange(2));
	trlstokeep=trldurs(trlstokeep,1);
	
	tempindex=ismember(msacstatsall2(:,1),trlstokeep);
	tempindex=find(tempindex);
	msacstatsall2=msacstatsall2(tempindex,:);
	
	tempindex=ismember(intstatsall2(:,1),trlstokeep);
	tempindex=find(tempindex);
	intstatsall2=intstatsall2(tempindex,:);
	
	rimsis2=intstatsall2(find(intstatsall2(:,3)==1),2);
	sacsts2=msacstatsall2(:,2);
	
else

	rimsis2=[];
	sacsts2=[];

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plotting Data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
	
axes(c(1))
text('Position',[1 1],'String',['saccade count:  ',num2str(msaccount)],'HorizontalAlignment','right')
text('Position',[1 .6],'String',['saccade rate:  ',num2str(msacrate),' sac/sec'],'HorizontalAlignment','right')
text('Position',[0 1],'String',['Saccade Statistics'],'FontSize',12)
text('Position',[0 .6],'String',['File:  ',savefilename],'FontSize',10)
text('Position',[0 .2],'String',['prog:  plotcompilemsac2.m'],'FontSize',8)
if ~isempty(rimsis2)
	text('Position',[1 .2],'String',['subset trl dur:  ',num2str(trldurrange(1)),' to ',num2str(trldurrange(2)),' msec'],'HorizontalAlignment','right')
end
axis('off')
	
axes(c(2))
maxbin=ceil(max(rimsis)/100)*100;
if ~isempty(xlimmaxes);maxbin=xlimmaxes(1);end
abins=[maxbin/50:maxbin/50:maxbin];
[n]=hist(rimsis,abins);
plot(abins,n/sum(n),'k-','LineWidth',2)
set(gca,'XLim',[0 maxbin*1.025],'XTick',[0:maxbin/10:maxbin])
if ~isempty(rimsis2)
	[n]=hist(rimsis2,abins);
	plot(abins,n/sum(n),'b--','LineWidth',1)
end
title('intersaccade intervals (real only; msec)','FontSize',10)

axes(c(3))
maxbin=ceil(max(sacsts)/10)*10;
if ~isempty(xlimmaxes);maxbin=xlimmaxes(2);end
abins=[maxbin/50:maxbin/50:maxbin];
[n]=hist(sacsts,abins);
plot(abins,n/sum(n),'k-','LineWidth',2)
set(gca,'XLim',[0 maxbin*1.025],'XTick',[0:maxbin/10:maxbin])
if ~isempty(sacsts2)
	[n]=hist(sacsts2,abins);
	plot(abins,n/sum(n),'b--','LineWidth',1)
end
title('saccade onset times (msec; from fp on)','FontSize',10)

