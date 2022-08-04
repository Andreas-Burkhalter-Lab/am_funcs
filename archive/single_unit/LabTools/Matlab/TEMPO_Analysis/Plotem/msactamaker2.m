function[raster,rasterin,rasterout]=msactamaker2(msac,emt,raster,filename,dispmode);

% function[msacta,msactatimes,msaccount]=msactamaker(msac,emt,raster,filename,dispmode);
% INPUTS:  msac is zeros and ones; emt is times and nans;
%	raster contains spike times in msec.
% CALLS:  rastermaker.m.

%% Setting Default Variables
eventnames={'All';
	    'Within Range';
	    'Outside Range'};

%%%%% Setting Root Configuration %%%%%
set(0,'ShowHiddenHandles','on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Checking Inputs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4;filename=[];end
if isempty(filename);filename=[];end
if nargin<5;dispmode=1;end
if isempty(dispmode);dispmode=1;end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figid='msactamaker.m';
uicolors=[.96 .57 .96;
	  .76 .54 1.0;
	  .65 .65 1.0;
	  .54 .76 1.0];

inclback=-300;
inclforw=300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Assessing Callback Status %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==1
   cbv=msac;
else
   cbv=0;
end

if cbv==0

	if dispmode
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Plotting Figure %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure('Units','inches','Position',[9 .3 5 6]);
  
	set(gcf,'Number','off','Name',figid)
	set(gcf,'DefaultAxesFontName','Palatino','DefaultAxesFontSize',8)
	set(gcf,'DefaultTextFontName','Palatino','DefaultTextFontSize',8)
	set(gcf,'DefaultAxesNextPlot','add','DefaultAxesUnits','inches')
	set(gcf,'DefaultAxesTickDir','out','DefaultUicontrolUnits','inches')
	set(gcf,'DefaultUicontrolFontName','Palatino','DefaultUicontrolFontSize',10)
	set(gcf,'PaperUnits','inches','PaperPosition',[1.75 2.25 8 10.5])

	axes('Position',[.5 .5 4 4])
 
    
    
   uicontrol('Style','pushbutton','Position',[3.7 5.5 .8 .28],'String', ....
      'CLOSE','Callback','msactamaker2(1);','FontWeight','bold')
   uicontrol('Style','pushbutton','Position',[2.8 5.5 .8 .28],'String', ....
      'SAVE','Callback','msactamaker2(2);','FontWeight','bold')
   s=[char(eventnames(1))];
   for x=2:size(eventnames,1)
		s=[s,'|',char(eventnames(x))];
   end   
   
	uic(1) = uicontrol('Style','popup','Position',[2.8 5.0 1.4 .28],'String',s, ....
		'BackGroundColor',uicolors(2,:).^.4)
	uic(2) = uicontrol('Style','edit','Position',[.5 5.0 0.6 .28], ....
		'HorizontalAlignment','left','FontSize',12, ....
    	'BackGroundColor',uicolors(1,:).^.4)
	uic(3) = uicontrol('Style','edit','Position',[1.8 5.0 0.6 .28], ....
      'HorizontalAlignment','left','FontSize',12, ...
      'BackGroundColor',uicolors(1,:).^.4)
   
   %% Headers
   text('Position',[-300 1.6],'String','Raster Output Parameters','FontSize',10);
   text('Position',[-300 1.45],'String','Start','FontSize',10);
	text('Position',[-100 1.45],'String','Range','FontSize',10);
	text('Position',[-200 1.30],'String','ms','FontSize',10);
	text('Position',[-8 1.30],'String','ms','FontSize',10);

	
     
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Getting Variables %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tottrl=size(emt,1);
	emtst=nanmin(emt')';
	emten=nanmax(emt')';
	msactatimes=[inclback:inclforw];
	msacta=zeros(1,length(msactatimes));
	msaccount=0;
   range = 25; %scan over 50 msec around each point during quantitative analysis  
 
   if ~isempty(raster)

		[x,y]=find(msac'==1);	% y = trl; x = index into emt column
		for sacnum=1:length(y)
			trl=y(sacnum);
			sacst=emt(trl,x(sacnum));
			if ((sacst+inclback)>=emtst(trl) & (sacst+inclforw)<=emten(trl))
				spikes=raster(trl,find(raster(trl,:)~=0))-sacst;
				index=find(spikes>=inclback & spikes<=inclforw);
				spikes=spikes(index);
				msacta([spikes-inclback+1])=msacta([spikes-inclback+1])+1;
				msaccount=msaccount+1;
			end
		end
	
		if msaccount>0
			msacta=msacta/msaccount;
		end
	
	end
   
   %msactatimes %time coordinates for each point
   %range of spike frequency for each point
   sacdist = [1:size(msacta)];  

	for limits = inclback:inclforw   
     	if (limits - range>inclback) & (limits + range<inclforw)
   		sacdist(find(msactatimes ==limits))=(sum(msacta(find(msactatimes == limits)-range:find(msactatimes == limits)+range))/(range*2+1));
      else sacdist(find(msactatimes ==limits)) = NaN;
      end
   end
   
   [sactrough, ttrough] = min(sacdist(find(msactatimes == 0): find(msactatimes == inclforw-2*range)));
   [saccrest, tcrest] = max(sacdist(find(msactatimes == 0): find(msactatimes == inclforw-2*range)));
   troughdist = msacta(find(msactatimes == ttrough)-range:find(msactatimes == ttrough)+range);
   crestdist = msacta(find(msactatimes == tcrest) - range:find(msactatimes == tcrest)+ range);
   basedist = msacta(find(msactatimes==(-2*range)):find(msactatimes ==(0)));
   
   
   rasterin =raster;
   rasterout=raster;
     
     
     
     
   %[sactrough, msactatimes(ttrough)]
   %[saccrest, msactatimes(tcrest)]
   meancrest = round(mean(crestdist)*1000);
   meantrough = round(mean(troughdist)*1000);
   meanbase = round(mean(basedist)*1000);
   
   if abs(meancrest - meanbase) > abs(meanbase - meantrough)
   	amplitude = round(meancrest - meanbase);   
   else amplitude = round(meantrough - meanbase);
   end
   
   disp('Significance between Max and Baseline')
   [h1,p1, ci1] = TTEST2(crestdist,basedist,.05,1); 
   %[p1, h1] = signrank(crestdist, basedist, .05)
   disp('Signficance between Min and Baseline')
   [h2,p2, ci2] = TTEST2(troughdist,basedist,.05,-1);
   %[p2, h2] = signrank(troughdist, basedist, .05)
   disp('Signficance between Max and Min')
   [h3,p3, ci3] = TTEST2(crestdist,troughdist,.05,1);
   %[p3, h3] = signrank(crestdist, troughdist,.05)
   
   %% 1st figure plot histograms and display data
   if dispmode
      plot(msactatimes,msacta*1000,'k-','LineWidth',.2)
		set(gca,'XLim',[min(msactatimes)-1 max(msactatimes)+1],'XTick',[min(msactatimes):50:max(msactatimes)])
		xlabel('msec');ylabel('spikes/second')
		ylims=get(gca,'YLim');yrange=ylims(2)-ylims(1);
		text('Position',[max(msactatimes) ylims(2)-.05*yrange],'String',['File:  ',filename],'FontSize',10,'HorizontalAlignment','right')
		title('Saccade-Triggered Spike Average','FontSize',12)

		y=normpdf([-6:6],0,4);y=y./(sum(y));
		z=conv2(msacta,y,'valid');
		offfront=ceil((length(y)-1)/2);
		offback=floor((length(y)-1)/2);
		z=[ones(1,offfront)*nan z ones(1,offback)*nan]; %z contains the spike/s averages.
      plot(msactatimes,sacdist*1000,'b-','LineWidth',2)
      plot(msactatimes,z*1000,'m-','LineWidth',2)	
	end
   
   %%% 2nd figure
   
      figure('Units','inches','Position',[9 .3 5 9])
		set(gcf,'Number','off','Name',figid)
		set(gcf,'DefaultAxesFontName','Palatino','DefaultAxesFontSize',8)
		set(gcf,'DefaultTextFontName','Palatino','DefaultTextFontSize',10)
		set(gcf,'DefaultAxesNextPlot','add','DefaultAxesUnits','inches')
		set(gcf,'DefaultAxesTickDir','out','DefaultUicontrolUnits','inches')
		set(gcf,'DefaultUicontrolFontName','Palatino','DefaultUicontrolFontSize',10)
		set(gcf,'PaperUnits','inches','PaperPosition',[1.75 2.25 8 10.5])

		axes('Position',[.5 7.8 4 1.5])
		axes('Position',[.5 5.5 4 2])
		axes('Position',[.5 3 4 2])
		axes('Position',[.5 .5 4 2])
      
      c=get(gcf,'children');
      c=sort(c);

  
		%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%% Plotting Data %%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%
	
		axes(c(7))
      text('Position',[0 0.6],'String',['Saccade Statistics    File:  ',filename],'FontSize',12)
      text('Position',[0 0.45],'String',['Minimum at ', num2str(ttrough),' ms:  Mean:', num2str(meantrough), ' spikes/s  Amplitude:', num2str(amplitude), ' spikes/s'])
      text('Position',[0 .3],'String',['Maximum at ', num2str(tcrest),' ms:  Mean:', num2str(meancrest), ' spikes/s'])
		text('Position',[0 .15],'String',['Baseline (-50 ms to 0 ms):  Mean:', num2str(meanbase), ' spikes/s'])
      text ('Position', [0 0], 'String', ['T-Test Stats:  Max and Baseline: ', num2str(p1), '  Min and Baseline: ' num2str(p2)])
      axis('off')
	
   	axes(c(8)) %for generating histograms of intersaccade intervals
  		
		%dum=ceil(max(troughdist(:,1))/10)*10
      hist (troughdist*1000,30)
		%set(gca,'XLim',[0 dum],'XTick',[0:10:dum])
		d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		title('Spike Frequency - Minimum','FontSize',10)

		axes(c(9))
		%dum=ceil(max(crestdist(:,4)))
		hist (crestdist*1000,30)
		%set(gca,'XLim',[0 dum],'XTick',[0:10:dum])
		d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		title('Spike Frequency - Maximum','FontSize',10)
      
		axes(c(10))
		%dum=ceil(max(msacstats(:,4)));
		%hist(msacstats(:,6),[.2:.2:dum])
      hist (basedist*1000,30)
		%set(gca,'XLim',[-180 180],'XTick',[-180:45:180])
		d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		title('Spike Frequency - Baseline','FontSize',10)

   
  
   
else

	if cbv == 1		%Close Button
		close(gcf)
	end
   
   if cbv == 2		%Save Raster File Button
      keyboard
   end
   
   
   
end