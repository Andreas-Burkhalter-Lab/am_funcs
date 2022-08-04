clear all;
[filename,path]=uigetfile('*.abf','pick your file');
cd (path); %set path
fullfilename=[path char(filename)];
data=abfload(fullfilename);
new=input('Is this from after 10/31/16?');
if ~new
    %Find start and stop of imaging
    inds=find((abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))))==1);
    meaninds=mean(diff(inds));
    figure;hold on;plot(data(:,5));plot(abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))),'r');
    figure;hold on;plot(data(:,5));plot(abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))),'r');
    xlim([inds(1)-2.5*meaninds inds(1)+2.5*meaninds]);
    [scanstart,y]=ginput(1);
    scanstart=round(scanstart);
    
    figure;hold on;plot(data(:,5));plot(abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))),'r');
    xlim([inds(end)-4*meaninds inds(end)+2*meaninds]);
    [scanstop,y]=ginput(1);
    scanstop=round(scanstop);
    disp(['Length of scan is ', num2str(scanstop-scanstart)])
    disp(['Time of scan is ', num2str((scanstop-scanstart)/1000)])
    
    close all;
    
    angle=data(scanstart:scanstop,4); %anglecut=view angle, channel 4 of "data"
    rewards=data(scanstart:scanstop,1);  %rewcut=rewards
    galvo=data(scanstart:scanstop,5);    %131018
    forwardvel=data(scanstart:scanstop,6);   %120806 EH
    rotationvel=data(scanstart:scanstop,7);  %120806 EH
    % plot(rewcut);ginput(1);
    ybinned=data(scanstart:scanstop,2);    %Ycut=yposition
    xbinned=data(scanstart:scanstop,3);
    numframes=length(scanstart:scanstop);
    if size(data,2)>7%131018
        ch8=data(scanstart:scanstop,8); %depending on pclamp protocol, may or may not have ch8
    end
else
    %Find start and stop of imaging
    inds=find((abs(diff(data(:,8)))>0.3*max(abs(diff(data(:,8)))))==1);
    meaninds=mean(diff(inds));
    figure;subplot(2,1,1);hold on;plot(data(:,8));plot(abs(diff(data(:,8)))>0.3*max(abs(diff(data(:,8)))),'r');
    subplot(2,1,2);hold on;plot(data(:,8));plot(abs(diff(data(:,8)))>0.3*max(abs(diff(data(:,8)))),'r');
%         xlim([scanstartrough-2.5*meaninds scanstartrough+2.5*meaninds]);
    xlim([inds(1)-2.5*meaninds inds(1)+2.5*meaninds]);
%         xlim([start.Position(1)-2.5 start.Position(1)+2.5])
    [scanstart,y]=ginput(1)
    scanstart=round(scanstart)
    
    plot(data(:,8));hold on;plot(abs(diff(data(:,8)))>0.3*max(abs(diff(data(:,8)))),'r');
%     xlim([scanstoprough-4*meaninds scanstoprough+2*meaninds]);
    xlim([inds(end)-4*meaninds inds(end)+2*meaninds]);
%         xlim([stop.Position(1)-2.5 stop.Position(1)+2.5])

    [scanstop,y]=ginput(1)
    scanstop=round(scanstop)
    disp(['Length of scan is ', num2str(scanstop-scanstart)])
    disp(['Time of scan is ', num2str((scanstop-scanstart)/1000/60),' mins'])
    
    close all;
    
    angle=data(scanstart:scanstop,4); %anglecut=view angle, channel 4 of "data"
    rewards=data(scanstart:scanstop,1);  %rewcut=rewards
    galvo=data(scanstart:scanstop,8);    %131018
    forwardvel=data(scanstart:scanstop,5);   %120806 EH
    rotationvel=data(scanstart:scanstop,6);  %120806 EH
    % plot(rewcut);ginput(1);
    ybinned=data(scanstart:scanstop,2);    %Ycut=yposition
    xbinned=data(scanstart:scanstop,3);
    numframes=length(scanstart:scanstop);
    
end

numfiles=input('Number of planes: ');
Ffile{numfiles}=0;
Ffilepath{numfiles}=0;

% % 170526, JS
rewards_th=1*rewards>4;
% % rewards=double(rewards_th);
rewards_df=diff(rewards_th);
rewards=[rewards_df(1) rewards_df']>=1;

% % 170526, commented by EBH because wrong
% rewards_th=1*rewards>(0.1*max(rewards));
% % % rewards=double(rewards_th);
% rewards_df=diff(rewards_th);
% rewards=[rewards_df(1) rewards_df']>=1;



%bin and average both forwardvel and rotationvel
%raw data need smoothing or binning to see well on compressed x scale.
%     velbinsize=200;   %  # of frames to bin.  50 looks ok.
%     binforwardvel=reshape(forwardvel,velbinsize,(numframes/velbinsize));    %gets 50 frames and puts into column.
%                                                                             %should have (numframes/50) columns
%     meanbinforwardvel=mean(binforwardvel);  %mean of each column(bin) and turns into vector of bins
%     binrotationvel=reshape(rotationvel,velbinsize,(numframes/velbinsize));    %for rotation
%     meanbinrotationvel=mean(binrotationvel);
%     timebinx=((velbinsize/2):velbinsize:(numframes-(velbinsize/2)));    %gives x(time) value of center of bin for plotting.
%
%     figure;
%     hold on;
%     plot(ybinned);
%     plot(rewards*600,'r');
%%
for n=1:numfiles
    [Ffile{n},Ffilepath{n}]=uigetfile('*.mat','pick the F file');
    fullFfile=[Ffilepath{n} Ffile{n}]
    %     load(fullFfile);
    %save(fullFfile,'ybinned','numframes','rewards','angle','forwardvel','rotationvel','velbinsize','meanbinforwardvel','meanbinrotationvel','timebinx','-append');
    pause(1);
    save(fullFfile,'ybinned','numframes','rewards','galvo','angle','forwardvel','rotationvel','xbinned','-append'); %131018 added galvobinned
    %     if size(data,2)>7%131018
    %         save(fullFfile,'ch8binned','-append');
    %     end
end

%%
% rew3=rewards+[rewards(2:end) 0]+[0 rewards(1:end-1)];
% rew=rew3(1:3:end);
endzonelength=10;
pos=(ybinned-min(ybinned));
pos=pos/max(pos)*180;
tops=pos>=(max(pos)-endzonelength);
bots=pos<=endzonelength;
topinds=find(tops);
botinds=find(bots);
% Easy way as long as there are only 2 out of place rewards and they are
%%not in end zones
rewlocs=find(rewards);
rewpossible=[topinds; botinds];
notExpected=rewlocs(~ismember(rewlocs,rewpossible));
% splitloc=find(notExpected);

% Harder way finding all rewards which look out of place
%
% tops = bwlabel(pos>(max(pos)-endzonelength));
% bots = -bwlabel(pos<endzonelength);
%
% % botcross=diff(bots);
% % botcross(botcross<0)=0;
% % botcross=-botcross;
% % topcross=diff(tops);
% % topcross(topcross<0)=0;
%
% crosses = tops+bots;
%
% inds = find(crosses);
% rewexpected = [];
% temp = -crosses(inds(1));
% for num = 1:length(inds)
%     if temp*crosses(inds(num))<0
%         rewexpected = [rewexpected crosses(inds(num))];
%         temp = crosses(inds(num));
%     end
% end
%
% notExpected = ~ismember(find(rew),rewexpected);

% Manually choose split points out of unexpected signals
disp(['Number of non-reward signals found: ',num2str(length(notExpected))]);
if isempty(notExpected)
    disp('No Environment Shifts Found')
    choice = input('Manually Add Points? [0-no/1-yes] ');
    if choice==1
%         if new==1
%            figure; plot(ybinned)
%            manstart=ginput(1);
%            manstart=manstart(1);
%            manstop=ginput(1);
%            manstop=manstop(1);
%         else
           figure; plot(ybinned)
           manstart=ginput(1);
           manstart=manstart(1);
           manstop=ginput(1);
           manstop=manstop(1);
%         end
%         manstart=input('Manually chosen start: ');
%         manstop=input('Manually chosen stop: ');
        timeNotExpected=[manstart manstop];
        timeSplit=timeNotExpected;
    end
    
%     for n=1:numfiles
%         fullFfile=[Ffilepath{n} Ffile{n}];
%         save(fullFfile,'timeSplit','-append'); %131018 added galvobinned
%     end
    
else
    timeOfSplits(length(notExpected))=0;
    timeNotExpected=(notExpected);
    figure; hold on;
    ignore = [];
    plot(pos)
    for i=1:length(timeNotExpected)
        line([timeNotExpected(i) timeNotExpected(i)], [min(pos) max(pos)],'Color','g','LineWidth',2)
    end
    for i=1:(max(length(timeNotExpected),2))
        if i<=length(timeNotExpected)
            line([timeNotExpected(i) timeNotExpected(i)], [min(pos) max(pos)],'Color','r','LineWidth',2)
        end
        choice = input('Accept split point [0-no/1-yes/2-add manual point/3-Manually set both points] ');
        if isempty(choice)|| choice == 0
            ignore=[ignore i];
        elseif choice==2 %This will only work for 2 split points. Need to deal with edge cases otherwise
            manrew=input('Manually chosen location: ');
            timeman=input('Before(0) or After(1) current?');
            if timeman==0
                timeNotExpected = [manrew timeNotExpected];
            else
                timeNotExpected = [timeNotExpected manrew];
            end
        elseif choice==3
            timeNotExpected=[];
           figure; plot(ybinned)
           manstart=ginput(1);
           manstart=manstart(1);
           manstop=ginput(1);
           manstop=manstop(1);
%         manstart=input('Manually chosen start: ');
%         manstop=input('Manually chosen stop: ');
        timeNotExpected=[manstart manstop];
        timeSplit=timeNotExpected;
            break;
        end
        
        
    end
    
    timeNotExpected(ignore)=[];
    
    % for split=1:length(notExpected)
    %     timeOfSplits(split)=(timeNotExpected(split));
    % end
    
    
    
    % split data into number of segments equal to number of non-reward linked
    %signals
    
    %  for i=1:(length(timeNotExpected)+1)
    %      if i==1
    %          start = scanstart;
    %          stop = timeNotExpected(i);
    %      elseif i == (length(timeNotExpected)+1)
    %          start = timeNotExpected(i-1)+1;
    %          stop = scanstop;
    %      else
    %          start = timeNotExpected(i-1)+1;
    %          stop = timeNotExpected(i);
    %      end
    % load(fullFfile);
    % vars=who('-file', fullFfile);
    % varst=cell2struct(vars,'var',length(vars));
    % newFfile=[fullFfile(1:(end-4)), 'Section',num2str(i),'.mat'];
    %
    %     angle=data(start:stop,4); %anglecut=view angle, channel 4 of "data"
    %     rewards=data(start:stop,1);  %rewcut=rewards
    %     galvo=data(start:stop,5);    %131018
    %     forwardvel=data(start:stop,6);   %120806 EH
    %     rotationvel=data(start:stop,7);  %120806 EH
    %    % plot(rewcut);ginput(1);
    %     ybinned=data(start:stop,2);    %Ycut=yposition
    %
    %     if size(data,2)>7%131018
    %         ch8=data(start:stop,8); %depending on pclamp protocol, may or may not have ch8
    %     end
    %
    %     ratio=length(F)/length(data);
    %     startF=ceil(ratio*start);
    %     stopF=ceil(ratio*stop);
    %     F=F(startF:stopF,:);
    %     Fc=Fc(startF:stopF,:);
    %     Fc2=Fc2(startF:stopF,:);
    %     Fc3=Fc3(startF:stopF,:);
    %     Fc3_DF=Fc3_DF(startF:stopF,:);
    %
    
    timeSplit=timeNotExpected;
    
    % ratio=round(length(F)/length(pos));
    % F=F(ratio*(start:stop),:);
    % F=F(ratio*(start:stop),:);
    % F=F(ratio*(start:stop),:);
    % F=F(ratio*(start:stop),:);
    % F=F(ratio*(start:stop),:);
    %
    for n=1:numfiles
        fullFfile=[Ffilepath{n} Ffile{n}];
        save(fullFfile,'timeSplit','-append'); %131018 added galvobinned
    end
    
end






