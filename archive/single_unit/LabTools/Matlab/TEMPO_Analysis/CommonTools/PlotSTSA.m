function[msacta,msactatimes,msaccount]=PlotSTSA(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol)    
% function[msacta,msactatimes,msaccount]=msactamaker(msac,emt,raster,filename,dispmode);
% INPUTS:  msac is zeros and ones; emt is times and nans;
%	raster contains spike times in msec.
% CALLS:  rastermaker.m.

Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

[msac, emt] = CalcMsac(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);    

raster = squeeze(data.spike_data(SpikeChan,:,BegTrial:EndTrial))';

%-----------------------------------------------------------------
%% TEMPORARY: SELECT SUBSET OF DATA FOR ZERO SPEED
%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );
%get the column of values of speeds in the dots_params matrix
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed(~null_trials)');
speed_select = logical( (speed > 0) );
msac = msac(speed_select, :);
emt = emt(speed_select, :);
raster = raster(speed_select, :);
if (length(unique_speed)>2) 
    disp(' more than two speeds not handled yet');
    return;
end
if (length(unique_speed)<2) 
    disp('This analysis only for runs with two speeds interleaved');
    return;
end
%-----------------------------------------------------------------

%%%%% Setting Root Configuration %%%%%
set(0,'ShowHiddenHandles','on')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figid='Saccade Triggered Spike Average';
uicolors=[.96 .57 .96;
    .76 .54 1.0;
    .65 .65 1.0;
    .54 .76 1.0];

inclback=-100;
inclforw=300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plotting Figure %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','inches','Position',[0 0 5 5.5]);

set(gcf,'Number','off','Name',figid)
set(gcf,'DefaultAxesFontName','Palatino','DefaultAxesFontSize',8)
set(gcf,'DefaultTextFontName','Palatino','DefaultTextFontSize',8)
set(gcf,'DefaultAxesNextPlot','add','DefaultAxesUnits','inches')
set(gcf,'DefaultAxesTickDir','out','DefaultUicontrolUnits','inches')
set(gcf,'DefaultUicontrolFontName','Palatino','DefaultUicontrolFontSize',10)
set(gcf,'PaperUnits','inches','PaperPosition',[1.75 2.25 8 10.5])

axes('Position',[.5 .5 4 4])



uicontrol('Style','pushbutton','Position',[3.7 4.8 .8 .28],'String', ....
    'CLOSE','Callback','msactamaker(1);','FontWeight','bold')	


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

[x,y]=find(msac'==1);	% y = trl; x = index into emt column
for sacnum=1:length(y)
    trl=y(sacnum);
    sacst(sacnum) = emt(trl,x(sacnum));
    if ((sacst(sacnum) +inclback)>=emtst(trl) & (sacst(sacnum) +inclforw)<=emten(trl))
        %count bins with multiple spikes accordingly
        spikes=sort([find(raster(trl,:) > 0 ) find(raster(trl,:) > 1 ) find(raster(trl,:) > 2 ) find(raster(trl,:) > 3 ) find(raster(trl,:) > 4 )]) - sacst(sacnum);
        index=find(spikes>=inclback & spikes<=inclforw);
        spikes=spikes(index); 
        
        %use these two lines to check
        sacstart(msaccount + 1) = sacst(sacnum);
        rast_data(msaccount+1,(spikes-inclback + 1)) = 1;
        for spike = 1: length(spikes)
            msacta([spikes(spike)-inclback+1])=msacta([spikes(spike)-inclback+1])+1;
        end   
        msaccount=msaccount+1;
    end
end

if msaccount>0
    msacta=msacta/msaccount;
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
basedist = msacta(find(msactatimes==(-range)):find(msactatimes ==(+range)));

meancrest = round(mean(crestdist)*1000);
meantrough = round(mean(troughdist)*1000);
meanbase = round(mean(basedist)*1000);

if abs(meancrest - meanbase) > abs(meanbase - meantrough)
    amplitude = round(meancrest - meanbase);   
else amplitude = round(meantrough - meanbase);
end

disp('Significance between Max and Baseline')
[h1,p1, ci1] = TTEST2(crestdist,basedist,.05,1); 
disp('Signficance between Min and Baseline')
[h2,p2, ci2] = TTEST2(troughdist,basedist,.05,-1);
disp('Signficance between Max and Min')
[h3,p3, ci3] = TTEST2(crestdist,troughdist,.05,1);

%% 1st figure plot 
plot(msactatimes,msacta*1000,'k-','LineWidth',.2)
set(gca,'XLim',[min(msactatimes)-1 max(msactatimes)+1],'XTick',[min(msactatimes):50:max(msactatimes)])
xlabel('msec');ylabel('spikes/second')
ylims=get(gca,'YLim');yrange=ylims(2)-ylims(1);
text('Position',[max(msactatimes) ylims(2)-.05*yrange],'String',['File:  ',FILE],'FontSize',10,'HorizontalAlignment','right')
title('Saccade-Triggered Spike Average','FontSize',12)

y=normpdf([-12:12],0,6);y=y./(sum(y));
z=conv2(msacta,y,'valid');
offfront=ceil((length(y)-1)/2);
offback=floor((length(y)-1)/2);
z=[ones(1,offfront)*nan z ones(1,offback)*nan]; %z contains the spike/s averages.
plot(msactatimes,sacdist*1000,'b-','LineWidth',2)
plot(msactatimes,z*1000,'m-','LineWidth',2)
legend('Raw', 'Boxcar Filter', 'Convolved', 2);
%boxcar filter with width of 2 * range


%print out summary data and append it to a log file
outfile = [BASE_PATH 'ProtocolSpecific\HDispTuning\STSA_summary.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE          base     crest    Pcrest   trough   Ptrough  ampl     ');
    fprintf(fid, '\r\n');
end
outstr = sprintf('%s %6.2f\t %6.2f\t %7.5f\t %6.2f\t %7.5f\t %6.2f\t ', FILE, meanbase, meancrest, p1, meantrough, p2, amplitude);
fprintf(fid, '%s', [outstr]);
fprintf(fid, '\r\n');
fclose(fid);

%append microsaccade times to  an array and store it to accumulate
fileid = [BASE_PATH 'ProtocolSpecific\HDispTuning\MicroSaccTimes.mat'];
if exist(fileid, 'file') ~= 0
    load (fileid)
else 
    num_files = 0;
    msac_times_out = [];
    sites = [];
end   

num_files = num_files + 1;
sites{num_files, 1} = FILE;      

sacstarts = find(msac == 1);
times = emt(sacstarts);
msac_times_out = [msac_times_out  times'];      

save (fileid, 'msac_times_out','num_files','sites');
return;