%% VR Part 1:  Import VR data
%Use cam (camera trigger signal) to break up other variables (stimdata
%corollary, 

%% Import Data and Organize

close all;
clear all;

inputmethod=input('Select by Dir(1) or File(2)?')

if inputmethod==1
    vrpath=uigetdir(pwd,'Pick VR Directory')
    cd(vrpath)
    allcsvfiles=ls('*.csv');
    for j=1:size(allcsvfiles,1)
        eval(['vrpath' num2str(j) ' = vrpath']);
        eval(['fullfilename' num2str(j) ' = [vrpath ''\'' allcsvfiles(j,:)]']);
        eval(['csvfilename' num2str(j) ' = allcsvfiles(j,:)']);
    end
    num_files=size(allcsvfiles,1);
    
elseif inputmethod==2
    num_files=input('Number of files?');

    for j=1:num_files
        [csvfilename,vrpath]=uigetfile('*.csv','pick your VR_*.csv file');
        cd (vrpath); %set path
        eval(['fullfilename' num2str(j) ' = [vrpath csvfilename]']);
         eval(['csvfilename' num2str(j) ' = csvfilename']);
    end
    
else
    error([inputmethod ' is an invalid input method']);    
end

keepers=who;
%%

for j=1:num_files
%j=1;
    cd(vrpath)
    eval(['fullfilename=fullfilename' num2str(j)]);
    eval(['csvfilename=csvfilename' num2str(j)]);
    
    sheet=1;
    [pathstr,name,ext]=fileparts(fullfilename);
    %[num,txt,raw] = xlsread(fullfilename,sheet);
    
    A=csvread(csvfilename,1,0);  %start on second row (1), exclude header; must be numeric
    time=A(:,1);  %time
    enc=A(:,2:3);  %two channels of rotary encoder
%             enc(:,1)=zeros(length(enc(:,1)),1);  %DELETE IF ENCODER CHANNEL IS BAD, ex. p12m1
    cam=A(:,4);  %camera recording start (rising edge) and stop (falling edge)
    stimdata=A(:,5);  %corollary copy of stimulus being presented
    
    savefn=strcat([name(4:end),'_vr.mat']);
%     save(savefn,'savefn','time','cam','enc','stimdata','name');
    

    %Determine trial begin and endpoints; parse each variable into trials
    clear trialstarts trialstops
    cam(2:size(cam,1),2)=diff(cam(:,1))>1;  %point-to-point difference reveals whether jump occurred
    cam(2:size(cam,1),3)=diff(cam(:,1))<-1;
    trialstarts(:,1)=find(cam(:,2)==1);  %jump up = frame number where trial (& camera recording) starts
    trialstops=find(cam(:,3)==1);  %jump down = frame number where trial (& camera recording) stops
    %Deal with occasional 'duplicate' count of a rise that is slower
    checkdupstarts=diff(trialstarts); checkdupstarts=[1000; checkdupstarts];
    checkdupstops=diff(trialstops); checkdupstops=[1000; checkdupstops];
    trialstarts(checkdupstarts<5)=[];
    trialstops(checkdupstops<5)=[];
    disp(length(trialstarts))
    disp(length(trialstops))


    
    if(trialstarts(1,1)>trialstops(1))
        error('Error: camera stopped before it started')
    elseif length(trialstops)==length(trialstarts)-1
        trialstops(end+1)=length(cam);
        disp('FYI, voltage recording was interrupted during final trial')
    end

    
    for t=1:length(trialstarts)
        stimdataT{t,1}=stimdata(trialstarts(t,1):min(trialstops(t,1),length(enc)));
        encT{t,1}=enc(trialstarts(t,1):min(trialstops(t,1),length(enc)),:);
        timeT{t,1}=time(trialstarts(t,1):min(trialstops(t,1),length(enc)),:);
    end
    
    temp=strsplit(name,'VR_');
    id=temp(2);
    
    id=name(4:end)
%% 
% confirm=input(['Validate id = ' id ' [Yes=1; No=Type id]'],'s');
% if confirm=='1'
%     id=id
% else
%     id=confirm
% end
%%
cd('Z:\2P Raw Data\BirthdatingFunctional\pCAGmCh3')
temp=strsplit(id,'m');
whichmap=strcat(temp{1},'stimMap');
load('Stimuli_13_3tri.mat',whichmap,'allStim');
eval(['[stim,stimVersion,flag_stimVersionNeedReconciled] = importStim(id,' whichmap ',allStim)']);
% eval(['stimVersion=' whichmap '(id)']);
% stimVersion=stimVersion(2);
stimOrig=stim;

flag_trialsNeedReconciled = (length(stimdataT)~=length(stim));
    
    cd(vrpath)
    save(savefn,'flag_stimVersionNeedReconciled','flag_trialsNeedReconciled','vrpath','id','stimVersion','stim','stimOrig','stimdataT','encT','timeT','time','A','trialstarts','trialstops','cam','enc','stimdata','name','savefn');
    
s=[];    
for i=1:length(keepers)
s=strcat(s,' ''',keepers{i},''' ')
end

eval(['clearvars -except ''keepers''' s]);
    
end