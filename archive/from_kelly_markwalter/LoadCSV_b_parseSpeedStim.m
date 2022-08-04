%% Load batch
clear;
close;

vrdir=uigetdir(pwd,'VR Directory') %('Z:\2P Raw Data\BirthdatingFunctional\pCAGmCh3\','VR Directory')
cd(vrdir)
allvrfiles=ls('*_vr.mat');
flag_someoneNeedsReconciling=0;


%%
for ka=6 %1:size(allvrfiles,1)
clearvars -except allvrfiles vrdir ka flag_someoneNeedsReconciling;

vrfile=allvrfiles(ka,:);
vrfileparts=strsplit(vrfile,'_vr.mat');
id=vrfileparts{1}

load(vrfile);   

%% Convert Encoder to RPM
    samprate=10; %kHz
    wind=10; %ms; time window to average speed
    
    %Determine timepoints of rising signal (1); flat or down signal (0)
    for t=1:length(trialstarts)
        encT{t}(2:size(encT{t},1),3)=diff(encT{t}(:,1))>0.1;  %Col3: Input 1 rising
        encT{t}(2:size(encT{t},1),4)=diff(encT{t}(:,2))>0.1;  %Col4: Input 2 rising
        encT{t}(2:size(encT{t},1),5)=diff(encT{t}(:,1))<-0.1;  %Col5: Input 1 falling
        encT{t}(2:size(encT{t},1),6)=diff(encT{t}(:,2))<-0.1;  %Col6: Input 2 falling
    end
    for t=1:length(trialstarts)
        encT{t}(2:size(encT{t},1),7)=diff(encT{t}(:,3))>0.5;  %Col 7: Input 1 rise begins (single pt)
        encT{t}(2:size(encT{t},1),8)=diff(encT{t}(:,4))>0.5;  %Col 8: Input 2 rise begins (single pt)
        encT{t}(2:size(encT{t},1),9)=diff(encT{t}(:,5))>0.5;  %Col 9: Input 1 fall begins (single pt)
        encT{t}(2:size(encT{t},1),10)=diff(encT{t}(:,6))>0.5;  %Col 10: Input 2 fall begins (single pt)
        encT{t}(:,11)=movsum(encT{t}(:,7),wind*samprate);  %Col 11: sum of counts in a moving window for # Input 1 rises
        encT{t}(:,12)=movsum(encT{t}(:,8),wind*samprate); %Col 12: sum of counts in a moving window for # Input 2 rises
        encT{t}(:,13)=movsum(encT{t}(:,9),wind*samprate); %Col 13: sum of counts in a moving window for # Input 1 falls
        encT{t}(:,14)=movsum(encT{t}(:,10),wind*samprate); %Col 14: sum of counts in a moving window for # Input 2 falls
    end 
        
    for t=1:length(trialstarts)
    speedsT{t}=999*ones(length(encT{t}),2); %make empty matrix (default value 999)
        for x = 1:length(encT{t})
           speedsT{t}(x,1)=wind*(x-1)/1000; %Time in ms(accounting for window)
           %Columns 11-15 of encT contain counts of number of transitions (up or down) for a channel (1 or 2); use average of the four to calculate speed
           %Units: 1000 ms/s; wind in ms; 60 s/min; 1024 steps/rev; num()= # steps in time window
           %Use this line if have only enc2 data, faulty enc1 data
%                 speedsT{t}(x,2)=(1000/wind)*(60/1024)*(encT{t}(x,12)+encT{t}(x,14))*0.5;  %averaging enc2 ups and downs; not using enc1, data bad
           %Use this line if have both enc1 and enc2 data
                speedsT{t}(x,2)=(1000/wind)*(60/1024)*(encT{t}(x,12)+encT{t}(x,14)+encT{t}(x,11)+encT{t}(x,13))*0.25;  %averaging enc1&2 ups and downs
        end  
    end
    
%% Add direction information to 'speed'; if not available, manually flip TreadRev trials
%Presuming sufficient encoder data, calculate whether direction is forward or reverse
    for t=1:length(trialstarts)
        direction=zeros(length(encT{t}),1);
        if max(sum(encT{t}(:,7:10),2))>1
            error('ERROR needs attention: Two encoder signals occur simultaneously at least once; Fwd vs Bkwd calculation may be affected');
        end

    direction(:,1)=encT{t}(:,7)*1;
    direction(:,2)=encT{t}(:,8)*2;
    direction(:,3)=encT{t}(:,9)*-1;
    direction(:,4)=encT{t}(:,10)*-2;
    direction(:,5)=sum(direction(:,1:4),2);

    [jrow,jcol,jumps]=find(direction(:,5)); %find non-zero elements --> row, col, value
    twoloc=find(jumps(1:length(jumps)-1)==2);  %jump# where enc2 rises
    twonext=jumps(twoloc+1);  %event immediately following enc2 rise event; should be 1 or -1 for fwd or bkwd, if -2, ignore and keep previous
    twotimes=timeT{t}(jrow(twoloc),1);

    irow=find(twonext~=-2);  idir=twonext(twonext~=-2); %eliminate mistakes when skipped

    direction(:,6)=zeros(length(direction(:,5)),1);
    for x=1:length(twonext)
    direction(jrow(twoloc(x)),6)=twonext(x);
    end
        
        % figure
        % subplot(2,1,1)
        % plot(timeT{t},speedsT{t}(:,2))
        % subplot(2,1,2)
        % plot(timeT{t},direction(:,6),'*')
        % ylim([-2 2])
        % suptitle(num2str(t))
        % ginput(1)
        % close

    direction(:,7)=direction(:,6);
    curr=1;
    direction(find(direction(:,6)==-2),6)=0;  %delete motor skips, ie when event -2 happens twice; will interpolate

    for x=1:length(direction(:,6))
        if direction(x,6)==0
            direction(x,7)=curr;
        else
            direction(x,7)=direction(x,6);
            curr=direction(x,7);
        end
    end

    speedsT{t}(:,3)=speedsT{t}(:,2).*direction(:,7);
    speedsT{t}(:,4)=smooth(speedsT{t}(:,3),500);  %Downsample 10kHz signal with 500 frame moving window --> 20Hz
    clear direction twoloc twonext twotimes irow idir curr jrow jcol jumps
    end
    save(savefn,'speedsT','-append');

%Check if either encoder channel is always zero, in which case directionality info is missing
 if max(enc(:,1)==0) || max(enc(:,2)==0)
    disp('WARNING: Lacking encoder channel data; no direction information');
    flag_NoEncDirectionInfo=true;
        %Then no directionality info (ex. missing an encoder); correct at least Bkwd Trials
        bkwdtrials=strfind(stim,'TreadReverse');
        bkwdtrials=find(~cellfun(@isempty,bkwdtrials));
    for idx=1:length(bkwdtrials)
        t=bkwdtrials(idx);
        speedsT{t}(20000:70000,4)=speedsT{t}(20000:70000,4)*-1;
    end
        %use stimDataT
%         figure
%         for idx=1:length(bkwdtrials)
%             subplot(3,1,idx)
%             t=bkwdtrials(idx)
%             plot(speedsT{t}(:,4),'r');
%             hold on
%             speedsT{t}(20000:70000,4)=speedsT{t}(20000:70000,4)*-1;
%             plot(speedsT{t}(:,4),'b');
%             suptitle('Bkwd Trials before (red) and after (blue) forced direction switch');
%             ylabel('Speed')
%         end

%Make concatenated version
speeds=[];
for k=1:length(speedsT)
    speeds=[speeds;speedsT{k}(:,4)];
end

        save(savefn,'speedsT','speeds','flag_NoEncDirectionInfo','-append');
 end


%%  % Figure - speed only, no stim data

% f=figure('units','normalized','outerposition',[0 0 1 1]);  %full-screen sized
% for t=1:length(trialstarts)
% subplot(6,7,t)
% plot(timeT{t},speedsT{t}(:,4),'r');
% title(['Trial #' num2str(t)])
% ylim([-30 30]);
% set(gca,'xtick',[]);
% xlim([timeT{t}(1) timeT{t}(end)])
% end
% suptitle([name ' - Speed by Trial (rpm)']);

%% Figure - speed and stimdata, all trials

%Figure for all trials
f=figure %('units','normalized','outerposition',[0 0 1 1]);  %full-screen sized
f.Position=[100 200 1100 700]
for t=1:length(stim)  %or alter to exclude BadAcq trial
subplot(6,7,t)
plot(timeT{t},stimdataT{t}*2-18,'r');
hold on
plot(timeT{t},speedsT{t}(:,4),'b');
text(timeT{t}(1)+1000,15,stim{t});
title(['Trial #' num2str(t)])
ylim([-20 20]);
set(gca,'xtick',[]);
xlim([timeT{t}(1) timeT{t}(end)])
end
namenounder=strrep(name,'_','-');
suptitle({[namenounder ' - Stimdata and Speeds as Recorded'], 'Speed Blue, Stimdata Red'});

cd(vrpath)
fn=[name '_speedsByTrial.jpg']
saveas(f,fn,'jpeg')
% fnfig=[name '_speedsByTrial.fig']
% saveas(f,fnfig)

%% Figure - Handle speed and stimdata, trial-by-trial to confirm

trialtypes=unique(stim);
g=figure
g.Position=[750 300 500 700]

% g=figure('units','normalized','outerposition',[0 0 1 1]);  %full-screen sized
for j=1:length(trialtypes)
    idxes=findTrialIndex(stim,trialtypes{j})
%     idxes=idxes(idxes~=39);  %to exclude BadAcq trial
    subplot(length(idxes)+1,1,1)
    ylabel('Stimulus');
    hold on
    for t=1:length(idxes)
        plot(timeT{idxes(t)}-timeT{idxes(t)}(1),stimdataT{idxes(t)}+t);
    end
    xlim([0 timeT{idxes(t)}(end)-timeT{idxes(t)}(1)])
for t=1:length(idxes)
    subplot(length(idxes)+1,1,t+1)
    plot(timeT{idxes(t)},stimdataT{idxes(t)}*2-18,'r');
    hold on
    plot(timeT{idxes(t)},speedsT{idxes(t)}(:,4),'b');
    ylim([-20 20]);
    set(gca,'xtick',[]);
    xlim([timeT{idxes(t)}(1) timeT{idxes(t)}(end)])
    title(['Trial #: ' num2str(idxes(t))]); 
end
suptitle(['Should be ' trialtypes{j} ' trials'])
fn=[name '_speedsByTrial_' trialtypes{j}];
saveas(g,fn,'jpeg');
ginput(1)
end

flag_TrialsNeedReconciled=input('TrialsNeedReconciled? [0=No, 1=Yes]');
cd(vrpath)
save(savefn,'flag_TrialsNeedReconciled','-append')
close(f);
close(g);

reconcileme='';
if flag_TrialsNeedReconciled==1
    flag_someoneNeedsReconciling=1;
    reconcileme=strcat(reconcileme,id,'&');
end

end

if flag_someoneNeedsReconciling==1
    disp(['WARNING: Need to Reconcile Trials of ' id])
else
    disp('Next Move on to IdentifyFramesToExclude.mat');
end

%%
%Saving
% fn=[name '_speedsByTrial.jpg']
% saveas(f,fn,'jpeg')
% fnfig=[name '_speedsByTrial.fig']
% saveas(f,fnfig)
% close

%save(savefn,'speedsT','-append');