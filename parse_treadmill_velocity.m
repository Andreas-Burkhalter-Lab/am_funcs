%% Load batch
% encT = two-column vector with inputs from the treadmill encoders recorded in Prairie View
% assumes that sample rate = 10khz
% %% Convert Encoder to RPM and meters per sec
%
%      speedsT = parse_treadmill_velocity(encT)
%
% updated 2019-04-17 on thermaltake

  
function speedsT = parse_treadmill_velocity(encT)

wheel_circumference_meters = 0.49; % circumference of treadmill used to convert revolutions per sec to meters per sec
samprate=10; %kHz
wind=10; %ms; time window to average speed
diff_threshold = 1; % threshold for detecting rises/falls in input channels; see notes to recording session 19048/2019-04-16 for details

%Determine timepoints of rising signal (1); flat or down signal (0)
encT(2:size(encT,1),3)=diff(encT(:,1)) > diff_threshold;  %Col3: Input 1 rising
encT(2:size(encT,1),4)=diff(encT(:,2)) > diff_threshold;  %Col4: Input 2 rising
encT(2:size(encT,1),5)=diff(encT(:,1)) < -diff_threshold;  %Col5: Input 1 falling
encT(2:size(encT,1),6)=diff(encT(:,2)) < -diff_threshold;  %Col6: Input 2 falling


encT(2:size(encT,1),7)=diff(encT(:,3))>0.5;  %Col 7: Input 1 rise begins (single pt) %%% acceleration (diff of diff) is greater than 0.5
encT(2:size(encT,1),8)=diff(encT(:,4))>0.5;  %Col 8: Input 2 rise begins (single pt)
encT(2:size(encT,1),9)=diff(encT(:,5))>0.5;  %Col 9: Input 1 fall begins (single pt)
encT(2:size(encT,1),10)=diff(encT(:,6))>0.5;  %Col 10: Input 2 fall begins (single pt)
encT(:,11)=movsum(encT(:,7),wind*samprate);  %Col 11: sum of counts in a moving window for # Input 1 rises
encT(:,12)=movsum(encT(:,8),wind*samprate); %Col 12: sum of counts in a moving window for # Input 2 rises
encT(:,13)=movsum(encT(:,9),wind*samprate); %Col 13: sum of counts in a moving window for # Input 1 falls
encT(:,14)=movsum(encT(:,10),wind*samprate); %Col 14: sum of counts in a moving window for # Input 2 falls



speedsT=999*ones(length(encT),2); %make empty matrix (default value 999)
    for x = 1:size(encT,1)
       speedsT(x,1)=wind*(x-1)/1000; %Time in ms(accounting for window) %%% 
       %Columns 11-15 of encT contain counts of number of transitions (up or down) for a channel (1 or 2); use average of the four to calculate speed
       %Units: 1000 ms/s; wind in ms; 60 s/min; 1024 steps/rev; num()= # steps in time window
       %Use this line if have only enc2 data, faulty enc1 data
%                 speedsT(x,2)=(1000/wind)*(60/1024)*(encT(x,12)+encT(x,14))*0.5;  %averaging enc2 ups and downs; not using enc1, data bad
       %Use this line if have both enc1 and enc2 data %%% multiply by 0.25 because one step = one up and one down each of Input 1 and Input 2
            speedsT(x,2)=(1000/wind)*(60/1024)*(encT(x,12)+encT(x,14)+encT(x,11)+encT(x,13))*0.25;  %averaging enc1&2 ups and downs %%% get revolutions per minute
    end  


%% Add direction information to 'speed'; if not available, manually flip TreadRev trials
%Presuming sufficient encoder data, calculate whether direction is forward or reverse

    direction=zeros(length(encT),1);
    if max(sum(encT(:,7:10),2))>1
        error('ERROR needs attention: Two encoder signals occur simultaneously at least once; Fwd vs Bkwd calculation may be affected');
    end

direction(:,1)=encT(:,7)*1;
direction(:,2)=encT(:,8)*2;
direction(:,3)=encT(:,9)*-1;
direction(:,4)=encT(:,10)*-2;    
direction(:,5)=sum(direction(:,1:4),2);

[jrow,jcol,jumps]=find(direction(:,5)); %find non-zero elements --> row, col, value
twoloc=find(jumps(1:length(jumps)-1)==2);  %jump# where enc2 rises
twonext=jumps(twoloc+1);  %event immediately following enc2 rise event; should be 1 or -1 for fwd or bkwd, if -2, ignore and keep previous
%     twotimes=timeT(jrow(twoloc),1);

irow=find(twonext~=-2);  idir=twonext(twonext~=-2); %eliminate mistakes when skipped

direction(:,6)=zeros(length(direction(:,5)),1);
for x=1:length(twonext)
direction(jrow(twoloc(x)),6)=twonext(x);
end

    % figure
    % subplot(2,1,1)
    % plot(timeT,speedsT(:,2))
    % subplot(2,1,2)
    % plot(timeT,direction(:,6),'*')
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

speedsT(:,3)=speedsT(:,2).*direction(:,7); % velocity
speedsT(:,4)=smooth(speedsT(:,3),500);  %Downsample 10kHz signal with 500 frame moving window --> 20Hz
speedsT(:,5) = speedsT(:,4) * wheel_circumference_meters / 60; % convert revolutions to meters, minutes to seconds (rpm to m/s)
clear direction twoloc twonext twotimes irow idir curr jrow jcol jumps

%     save(savefn,'speedsT','-append');

%% commented out following section - causes error when wheel data was not recorded 
% %Check if either encoder channel is always zero, in which case directionality info is missing
%  if max(encT(:,1)==0) || max(encT(:,2)==0)
%     disp('WARNING: Lacking encoder channel data; no direction information');
%     flag_NoEncDirectionInfo=true;
%         %Then no directionality info (ex. missing an encoder); correct at least Bkwd Trials
%         bkwdtrials=strfind(stim,'TreadReverse');
%         bkwdtrials=find(~cellfun(@isempty,bkwdtrials));
%     for idx=1:length(bkwdtrials)
%         t=bkwdtrials(idx);
%         speedsT(20000:70000,4)=speedsT(20000:70000,4)*-1;
%     end
%     
%     
% 
%     
%         %use stimDataT
% %         figure
% %         for idx=1:length(bkwdtrials)
% %             subplot(3,1,idx)
% %             t=bkwdtrials(idx)
% %             plot(speedsT(:,4),'r');
% %             hold on
% %             speedsT(20000:70000,4)=speedsT(20000:70000,4)*-1;
% %             plot(speedsT(:,4),'b');
% %             suptitle('Bkwd Trials before (red) and after (blue) forced direction switch');
% %             ylabel('Speed')
% %         end
% 
% %Make concatenated version
% % speeds=[];
% % for k=1:size(speedsT,1)
% %     speeds=[speeds;speedsT{k}(:,4)];
% % end
% 
% %         save(savefn,'speedsT','speeds','flag_NoEncDirectionInfo','-append');
%  end
%%
 
 end
%%  % Figure - speed only, no stim data

% f=figure('units','normalized','outerposition',[0 0 1 1]);  %full-screen sized
% 
%     subplot(6,7,t)
%     plot(timeT,speedsT(:,4),'r');
%     title(['Trial #' num2str(t)])
%     ylim([-30 30]);
%     set(gca,'xtick',[]);
%     xlim([timeT(1) timeT(end)])
% 
% suptitle([name ' - Speed by Trial (rpm)']);