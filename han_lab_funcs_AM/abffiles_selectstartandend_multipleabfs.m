%131018 added galvo and lick sensor import. not tested yet

clear all;
close all;

%abf_files: ch1=rewards; ch2=y-position; ch3=x-position; ch4=view angle; ch5=slow galvo; ch6=forward ball velocity; ch7=rotation ball velocity

% numfiles=3;
% numframes=24000; %total for all files

%added input dialog to get numfiles and frames.  120815  EH
prompts={'Number of .abf files','Number of total frames'};
ptitle=('File info');
lineno=1;
defaultanswers={'11','28000'};
inp=inputdlg(prompts, ptitle, lineno,defaultanswers);
numfiles=str2num(inp{1});
numframes=str2num(inp{2});


    
ybinned=[];
rewards=[];
angle=[]; %for view angle.  EH 120806
forwardvel=[];  %forward velocity.  EH 120806
rotationvel=[]; %rotational velocity.  EH 120806
ch8binned=[];   %131018
galvobinned=[]; %131018
for jj=1:numfiles
    ybinnedtemp=[];
    rewtemp=[];
        galvobinnedtemp=[];%131018
    ch8binnedtemp=[];%131018
    [filename,path]=uigetfile('*.abf','pick your files');
    cd (path); %set path
                %for moving all scripts to "Local"  120628  EH
    fullfilename=[path char(filename)];
    data=abfload(fullfilename);
    
    
   % data=data(1.261e5:3.248e5,:); %cut out part of abf trace
    
    inds=find((abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))))==1);
      meaninds=mean(diff(inds));
      figure;hold on;plot(data(:,5));plot(abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))),'r');
     figure;hold on;plot(data(:,5));plot(abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))),'r');
    xlim([inds(1)-2.5*meaninds inds(1)+2.5*meaninds]);
    [scanstart,y]=ginput(1)
    scanstart=round(scanstart)
    
     figure;hold on;plot(data(:,5));plot(abs(diff(data(:,5)))>0.3*max(abs(diff(data(:,5)))),'r');
    xlim([inds(end)-4*meaninds inds(end)+2*meaninds]);
    [scanstop,y]=ginput(1)
    scanstop=round(scanstop)
    disp(['Length of scan is ', num2str(scanstop-scanstart)])
%     if jj==1
%        scanstart=15113;
%         scanstop=445184;
%     end
%     

%         if jj==2
%         scanstart=2583;
%         scanstop=217620;
%         end


%         if jj==3
%         scanstart= 6814;
%          scanstop=221851;
%         end
% 
%         if jj==4
%         scanstart=7356;
%          scanstop=222391;
%         end
% 
%         if jj==5
%         scanstart=13737;
%          scanstop=228774;
%         end
% 
%         if jj==6
%         scanstart=5442;
%          scanstop=220479;
%         end
   
    close all;

    anglecut=data(scanstart:scanstop,4); %anglecut=view angle, channel 4 of "data"
    rewcut=data(scanstart:scanstop,1);  %rewcut=rewards
    galvocut=data(scanstart:scanstop,5);    %131018
    forwardvelcut=data(scanstart:scanstop,6);   %120806 EH
    rotationvelcut=data(scanstart:scanstop,7);  %120806 EH
   % plot(rewcut);ginput(1);
    Ycut=data(scanstart:scanstop,2);    %Ycut=yposition
    
    if size(data,2)>7%131018
        ch8cut=data(scanstart:scanstop,8); %depending on pclamp protocol, may or may not have ch8
    end
    
    binsize=numfiles*length(Ycut)/(numframes);
    %bin Y and angle 
    for i=1:numframes/numfiles
        ybinnedtemp(i)=mean(Ycut(round((i-1)*binsize+1):floor(round(i*binsize))));
        rewtemp(i)=mean(rewcut(round((i-1)*binsize+1):floor(round(i*binsize))));
        angletemp(i)=mean(anglecut(round((i-1)*binsize+1):floor(round(i*binsize)))); %120806 EH
        forwardveltemp(i)=mean(forwardvelcut(round((i-1)*binsize+1):floor(round(i*binsize)))); %120806 EH
        rotationveltemp(i)=mean(rotationvelcut(round((i-1)*binsize+1):floor(round(i*binsize)))); %120806 EH
        galvobinnedtemp(i)=mean(galvocut(round((i-1)*binsize+1):floor(round(i*binsize))));%131018
            if size(data,2)>7%131018
                ch8binnedtemp(i)=mean(ch8cut(round((i-1)*binsize+1):floor(round(i*binsize))));
            end
    end
    ybinned=[ybinned ybinnedtemp];
    rewards=[rewards rewtemp];
    angle=[angle angletemp];    %120806 EH
    forwardvel=[forwardvel forwardveltemp]; %120806 EH
    rotationvel=[rotationvel rotationveltemp];  %120806 EH
    galvobinned=[galvobinned galvobinnedtemp];  %131018
    
        if size(data,2)>7%131018
            ch8binned=[ch8binned ch8binnedtemp];
        end
    
    
end

     rewards_th=rewards>(0.1*max(rewards));
    rewards_df=diff(rewards_th);
    rewards=[rewards_df(1) rewards_df]>=1;
    
    %bin and average both forwardvel and rotationvel
    %raw data need smoothing or binning to see well on compressed x scale.
%     velbinsize=200;   %  # of frames to bin.  50 looks ok.
    velbinsize=gcd(numframes,200); %uhhhh mwa
    binforwardvel=reshape(forwardvel,velbinsize,(numframes/velbinsize));    %gets 50 frames and puts into column. 
                                                                            %should have (numframes/50) columns
    meanbinforwardvel=mean(binforwardvel);  %mean of each column(bin) and turns into vector of bins
    binrotationvel=reshape(rotationvel,velbinsize,(numframes/velbinsize));    %for rotation
    meanbinrotationvel=mean(binrotationvel);
    timebinx=((velbinsize/2):velbinsize:(numframes-(velbinsize/2)));    %gives x(time) value of center of bin for plotting.
 
    figure;
    hold on;
    plot(ybinned);
    plot(rewards*600,'r');
for p=1:4
[Ffile,Ffilepath]=uigetfile('*.mat','pick the F file');
fullFfile=[Ffilepath Ffile];

 %save(fullFfile,'ybinned','numframes','rewards','angle','forwardvel','rotationvel','velbinsize','meanbinforwardvel','meanbinrotationvel','timebinx','-append'); 
save(fullFfile,'ybinned','numframes','rewards','galvobinned','angle','forwardvel','rotationvel','-append'); %131018 added galvobinned
 if size(data,2)>7%131018
     save(fullFfile,'ch8binned','-append');
 end
end