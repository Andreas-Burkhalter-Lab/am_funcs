function[summarystats]=plotmsacstats(msac,emt,emh,emv,index,datapath,filename,theta,trlnums,raster1, foo,dispmode);

% function[summarystats]=plotmsacstats(msac,emt,emh,emv,index,filename,theta,dispmode)
% Function to plot microsaccade statistics.
% INPUTS:  msac is zeros and ones; emt is times and nans; emh and emv are
%	eye position values in degrees; filename is optional and is placed
%	at the top of the figure if supplied.
%	trlnums = absolute trl numbers of all data
% OUTPUT:  summarystats contains --
%	[msaccount msacrate msacdurmean msacdurstd msacamplrmean msacamplrstd]

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Checking Inputs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<6;filename=[];end
if isempty(filename);filename=[];end
if nargin<12;dispmode=1;end
if isempty(dispmode);dispmode=1;end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figid='plotmsacstats.m';
uicolors=[.96 .57 .96;
	  .76 .54 1.0;
	  .65 .65 1.0;
	  .54 .76 1.0];

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
		figure('Units','inches','Position',[9 .3 5 9.5])
		set(gcf,'Number','off','Name',figid)
		set(gcf,'DefaultAxesFontName','Palatino','DefaultAxesFontSize',8)
		set(gcf,'DefaultTextFontName','Palatino','DefaultTextFontSize',10)
		set(gcf,'DefaultAxesNextPlot','add','DefaultAxesUnits','inches')
		set(gcf,'DefaultAxesTickDir','out','DefaultUicontrolUnits','inches')
		set(gcf,'DefaultUicontrolFontName','Palatino','DefaultUicontrolFontSize',10)
		set(gcf,'PaperUnits','inches','PaperPosition',[1.75 2.25 8 10.5])

		axes('Position',[.5 8.3 4 .5])
		axes('Position',[.5 6.5 4 1.5])
      axes('Position',[.5 4.5 4 1.5])
      axes('Position',[.5 2.5 4 1.5])
      axes('Position',[.5 .5 4 1.5])
     		uicontrol('Style','pushbutton','Position',[3.7 9 .8 .28],'String', ....
			'CLOSE','Callback','plotmsacstats(1);','FontWeight','bold')
      c=get(gcf,'children');
      c=sort(c);
	
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Getting Variables %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	[exptstats,trlstats,msacstats,intstats]=getmsacstats(msac,emt,emh,emv);
	msaccount=exptstats(1);
	msacrate=exptstats(2);
	msacdurmean=mean(msacstats(:,3));
   msacdurstd=std(msacstats(:,3));
   msacdir=(msacstats(:,6));
   
   for sac = 1: size(msacstats,1)
      sacstart(sac) = msacstats(sac,1) - index(msacstats(sac,5), 3);
   end 
   
   
   
   
   if theta < 90 | theta > 270 	% calculate difference in direction based on system of -180 to +180
      if theta > 180
         theta = theta - 360;		% convert theta to system from -180 to +180
      end
      thetadiff = abs(mean(msacdir) - theta);
      
      for dummy=1:size(msacdir)	%%Change from values or -180 -> 180 to 0 -> 360
   		if msacdir(dummy)<0
         	msacdir(dummy)=msacdir(dummy)+360;
         end   
      end
   else % calculate difference in direction based on system from 0 to 360
      for dummy=1:size(msacdir)	%%Change from values or -180 -> 180 to 0 -> 360
      	if msacdir(dummy)<0
         	msacdir(dummy)=msacdir(dummy)+360;
      	end
      end
      thetadiff = abs(mean(msacdir) - theta);
   end
   
   msacampmean=mean(msacstats(:,4));
	msacampstd=std(msacstats(:,4));
	rimsis=intstats(find(intstats(:,3)==1),1);
   pimsis=intstats(find(intstats(:,3)==0),1);
   
     
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Creating Output %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	summarystats=[msaccount msacrate msacdurmean msacdurstd msacampmean msacampstd];
	
	if dispmode
	
		%%% Displaying Results
		disp(' ')
		disp(['saccade rate:  ',num2str(msacrate),' saccades/second'])
		disp(['saccade duration:  ',num2str(msacdurmean),' (+- ',num2str(msacdurstd),') msec'])
		disp(['saccade amplitude:  ',num2str(msacampmean),' (+- ',num2str(msacampstd),') degrees'])
		disp(' ')
	
		%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%% Plotting Data %%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%
	
		axes(c(6))
		text('Position',[1 1],'String',['saccade count:  ',num2str(msaccount)],'HorizontalAlignment','right')
		text('Position',[1 .6],'String',['saccade rate:  ',num2str(msacrate),' sac/sec'],'HorizontalAlignment','right')
		text('Position',[0 1],'String',['Saccade Statistics'],'FontSize',12)
		text('Position',[0 .6],'String',['File:  ',filename],'FontSize',10)
      text ('Position', [0 0], 'String', ['Median for intersaccade intervals= ', num2str(median(intstats(:,1))), ....
            ' Mean Dir - Pref Dir: ', num2str(thetadiff)], 'HorizontalAlignment','left')
      axis('off')
	
   	axes(c(7)) %for generating histograms of saccade start time - stimulus onset time
      dum=ceil(max(sacstart)/10)*10;
      hist(sacstart,[-1000:50:dum])
      set(gca,'XLim',[0 dum],'XTick',[0:100:dum])
		d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		title('Time of Saccade offset following Stimulus Onset (msec)','FontSize',10)
      
      %axes(c(7)) %for generating histograms of intersaccade intervals
    	%dum=ceil(max(intstats(:,1))/10)*10;
      %hist(intstats(:,1),[0:50:dum])
		%set(gca,'XLim',[0 dum],'XTick',[0:100:dum])
		%d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		%title('intersaccade intervals (msec)','FontSize',10)
      
      
      
      %axes(c(7)) %for generating histograms of saccade durations
		%dum=ceil(max(msacstats(:,3))/10)*10;
		%hist(msacstats(:,3),[5:5:dum])
		%set(gca,'XLim',[0 dum],'XTick',[0:5:dum])
		%d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		%title('saccade durations (msec)','FontSize',10)
      
		axes(c(8))
		dum=ceil(max(msacstats(:,4)));
		hist(msacstats(:,4),[.2:.2:dum])
		set(gca,'XLim',[0 dum],'XTick',[0:.4:dum])
		d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		title('saccade vector amplitudes (degrees)','FontSize',10)
      
		axes(c(9))
		%dum=ceil(max(msacstats(:,4)));
		%hist(msacstats(:,6),[.2:.2:dum])
      hist(msacdir,20)
		set(gca,'XLim',[0 360],'XTick',[0:45:360])
		d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		title('saccade directions (degrees)','FontSize',10)
      
      axes(c(10))
		%dum=ceil(max(msacstats(:,4)));
		%hist(msacstats(:,6),[.2:.2:dum])
      plot(msacstats(:,4), msacstats(:,7), 'o')
		%set(gca,'XLim',[0 360],'XTick',[0:45:360])
	%	d=get(gca,'children');set(d(5),'FaceColor',[.7 .7 .7]);
		title('Main Sequence Analysis','FontSize',10)

      
      
	end
   
   output = 0;
   
   if output == 1
   %save saccade stats in file
	  if size(datapath,1) == 1 & size(datapath,2) == 1		%convert PATH to string if PATH is cell array
	  		datapath = char(datapath);		
	 end
  	 if size(filename,1) == 1 & size(filename,2) == 1		%convert PATH to string if PATH is cell array
  			filename = char(filename);		
  	 end
   
 	  i = size(datapath,2) - 1;
	   while datapath(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
 	     i = i - 1;
	   end   
  	 datapathout = [datapath(1:i) 'Analysis\Saccades\'];
   
	   %add an extension to the filename
  	 i = size(filename,2);
  	 fileout = [filename(1:i) '.sac'];
   
   
		rastoutstart = -100;	%range in ms
  	 	rastoutend = +300;
   
   	disp ('Outputting spike times');   
	   fileid = [datapathout fileout];
 		fwriteid = eval(['fopen(fileid, ''w'')']);
   	for sac = 1: size(msacstats,1)
   	   %difference between this and msactamaker - key off of end of saccade here if use msacstats (sac,2)
	   	if ( msacstats(sac,1) + rastoutstart >= index(msacstats(sac,5), 3) ) & ( msacstats(sac,1) + rastoutend <= index(msacstats(sac,5), 4) )
    	     fprintf(fwriteid, '%i	%.2f	%.2f	%i	%i	%.4f	%6.2f	%.2f\n', trlnums(msacstats(sac,5)), foo(msacstats(sac,5), 1), foo(msacstats(sac,5),2), msacstats(sac,1), msacstats(sac,3), msacstats(sac,4), msacstats(sac,7), msacdir(sac));	
    	  end   
  		end    
 	   fclose(fwriteid);
   
 	  %I added the next set of code to expedite the analysis and generation of PSTAs.  - BJP 3/21/00
	  %output peri-saccade spike times
    
   i = size(datapath,2) - 1;
   while datapath(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
      i = i - 1;
   end   
   datapathout = [datapath(1:i) 'Analysis\STSA\'];
   
   disp ('Outputting Saccade Stats');
   %add an extension to the filename
   i = size(filename,2);
   fileout = [filename(1:i) '.sta'];
     
   fileid = [datapathout fileout];
   fwriteid = eval(['fopen(fileid, ''w'')']);
      
   for sac = 1: size(msacstats,1)
     	%difference between this and msactamaker - key off of end of saccade here if use msacstats (sac,2)
	   if ( msacstats(sac,1) + rastoutstart >= index(msacstats(sac,5), 3) ) & ( msacstats(sac,1) + rastoutend <= index(msacstats(sac,5), 4) )
      	fprintf(fwriteid, '%i	%i	', trlnums(msacstats(sac,5)), msacstats(sac,1));	
      	sacstart = msacstats(sac,1);
      
      	%section out appropriate spike times from raster
      	spikevect = raster1(msacstats(sac,5),:);
      	spikevect = spikevect((sacstart + rastoutstart) <= spikevect & spikevect <= (sacstart + rastoutend));
      	% next line offsets spike times to be time 0 at start of saccade   
      
      	%these lines are for checking output of spike times
      	%rastout = (spikevect - sacstart) - rastoutstart +1;
      	%rast_dat(sac, rastout) = 1;
      
      	rastout = spikevect - sacstart;
      
      	fprintf(fwriteid,'%i ',rastout);
      	fprintf(fwriteid,'\n');
      end
	end    
   fclose(fwriteid);

	end %output = 0 or 1   
   
else

	if cbv == 1
		close(gcf)
	end

end
