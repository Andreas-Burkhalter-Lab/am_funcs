%To Do:
%In "all_events", find bin max from frame max
%
%In making ave events, skips traces with NaNs (cut off traces).  Want to
%get those traces into averages.

%Place vs. non place events.

clear all;
close all;
num_cells=[];
num_fields=[];
field_sizes=[];
fs=31;
planes=6;
ms_per_frame=1000 / (fs/planes);
frame_rate_corr= 742.4 / ms_per_frame;
% frame_rate_corr=30;
%ms_per_frame=71.68; %71.68 for 128 lines X 0.56ms
%18.56 for 32 lines x 0.56
%frame_rate_corr=10; %40 for 18.56 ms /frame images
%20 for 37.12 ms /frame images
%10 for 74.24 ms /frame images
%used in finding pos or neg long runs.
%allows x frames to be in opposite direction without
%throwing out entire run.
%number has to be higher at higher frame rate
%based on empirical testing so feel free to change


%120815 EH
% h=warndlg('Did you enter proper F-file name in script?');%reminder to paste in F-file below
% waitfor(h);

for iii=1:1
    switch iii
        case 1
            %            load 120817_054_motcorr_roibyhand_F.mat; zoom=7; Fname='120817_054_motcorr_roibyhand_F.mat'; %17,18,19,20
            [Ffile,Ffilepath]=uigetfile('*.mat','pick the F file');
            [pathstr, name, ext] = fileparts(Ffile);    %need to get stem of name without ext for saving figs
            cd (Ffilepath); %set path
            %for moving all scripts to "Local"  120628  EH
            fullFfile=[Ffilepath Ffile];
            load(fullFfile);
            
            
            
    end
    
    ms_per_line=0.5;
    
    minfieldwidth=90;   %120 originally, i think.
    minratio=3; %3 originally
    thresh1=0.25;
    minDF=0.1;
    minrate=0.3; %0.3 originally
    baselength=0.25;
    
    combinedYpos_rew=[];
    combinedYpos_no_rew=[];
    
    combinedYneg_rew=[];
    combinedYneg_no_rew=[];
    
    
    vthresh=2;%4 originally
    ysize_thresh=200; %400 originally
    
    trackstart=-850;
    trackend=850;
    tracksize=trackend-trackstart;
    
    
    numneurons=size(Fc3,2);
    
    [~,ybinned,~,~,~,~,~,rewards,~]= align_running({Ffilepath},{Ffile},1);
    %%
    ybinned=ybinned(:,1);
        [~,reward_ind]=findpeaks(rewards,'MinPeakDistance',4);   %rewards (logic) in frame# where 1=reward
    figure;

    hold on;
    plot(ybinned);  %ybinned is Ypos.  binning?
%     plot(rewards*1600-800);
for rr=1:length(reward_ind)
   line([(reward_ind(rr)) (reward_ind(rr))],[min(ybinned) max(ybinned)],'color','r')
end
    
    
    ydiff=diff(ybinned);    %difference between adjacent ypos.  speed.
    ydiff=[ydiff(1); ydiff]; %duplicates 1st pos. so will have same # of samples as before.
    
    %positive direction
    L = bwlabel(ydiff>vthresh); %vector where each run faster  than
    %vthresh is marked.
    %successive events with incrementing #s,
    %i.e. 1,2,3, etc.
    %"L" gets reused for different things
    ysize=[];
    
    
    
    for i=1:max(L)  %gets ysize of each run
        ysize(i)=abs(ybinned(find(L==i,1,'first'))-ybinned(find(L==i,1,'last')));
        if ysize(i)<ysize_thresh
            L(L==i)=0;  %turns events shorter than ysize_thresh to 0.
        end
    end
    L_pos=L>0;  %L has ascending #s for contig runs from bwlabel
    %L_pos is logical with 1 at frames where longrun, pos.
    longruns_pos=find(L>0); %longruns_pos=frame #s where pos, long run
    
    %negative direction
    L = bwlabel(ydiff<-vthresh);
    ysize=[];
    for i=1:max(L)
        ysize(i)=abs(ybinned(find(L==i,1,'last'))-ybinned(find(L==i,1,'first')));
        if ysize(i)<ysize_thresh
            L(L==i)=0;
        end
    end
    L_neg=L>0;
    longruns_neg=find(L>0);
    
    %for bringing in Dan's placecells_longruns script.  121024  EBH
    %gets peaks and widths of events during any long run
    %second part is in nested "if" statements.
    %     longrunstrace=union(longruns_pos,longruns_neg);
    %     longruns=find(L>0);
    %modified because want only longruns with reward
    
    
    %turn Fc3 back into DF/F (from std)
    Fc3_DF=zeros(size(Fc2));
    event_dur=[];
    event_summary_all=[];
    Fc3_DF_events_label=[];
    duration=(length(ybinned)*ms_per_frame)/1000;   %duration of imaging (s)
    for i=1:numneurons
        upperbase(i)=3*std(Fc2(:,i));
        baselineind=find(Fc2(:,i)<upperbase(i));
        baseline=Fc2(baselineind,i);
        basemedian=median(baseline);
        basestd=std(baseline);
        Fc3_DF(:,i)=Fc3(:,i)*basestd;
        %             %get information on Fc3_DF events.  120919 EH
        %
        %             Fc3_DF_events_label(i)=bwlabel(Fc3_DF(:,i),4);
        %             num_Fc3_events=max(Fc3_DF_events_label(i));
        %             event_freq=num_Fc3_events/duration;
        %             %event_dur=zeros(i,num_Fc3_events);
        %
        %                 for events_loop=1:num_Fc3_events
        %                     event_dur(events_loop)=(length(Fc3_DF_events_label(find(Fc3_DF_events_label==events_loop,i))))*(ms_per_frame/1000);
        %                 end
        %               %one row for each neuron,...
        %               %columns: # of events, event freq (s), event dur (ms)
        %             event_summary_all(i,1)= num_Fc3_events;
        %             event_summary_all(i,2)= event_freq;
        %             event_summary_all(i,3)= mean(event_dur);
        %             event_summary_all(i,4)= std(event_dur);
        %             clear Fc3_DF_events_label num_Fc3_events event_freq  event_dur;
        
    end
    
    %event_summary_all%prints info
    
    %break long runs into reward or not
    %first find reward runs
    
    
  
    
    %now find all long runs that are not reward runs
    longruns_neg_no_rew=[];
    for tt=1:length(longruns_neg)
        if sum(find(longruns_neg_rew==longruns_neg(tt)))==0
            longruns_neg_no_rew=[longruns_neg_no_rew; longruns_neg(tt)];
        end
        
    end
    
    longruns_pos_no_rew=[];
    for tt=1:length(longruns_pos)
        if sum(find(longruns_pos_rew==longruns_pos(tt)))==0
            longruns_pos_no_rew=[longruns_pos_no_rew; longruns_pos(tt)];
        end
        
    end
    
    
    

    test_neg_no=zeros(length(Fc3),1);
    for tt=1:length(longruns_neg_no_rew)
        test_neg_no(longruns_neg_no_rew(tt))=1;
    end
    
    test_pos_no=zeros(length(Fc3),1);
    for tt=1:length(longruns_pos_no_rew)
        test_pos_no(longruns_pos_no_rew(tt))=1;
    end
    
    ybinnedt=ybinned-min(ybinned);
    figure; title ('pos'); hold on; plot(ybinnedt/max(ybinnedt),'r');plot(test_pos_no,'g');
    
    figure; title ('neg'); hold on; plot(ybinnedt/max(ybinnedt),'r');plot(test_neg_no,'g');
    
    
    
    %bin position and average Ftrace
    
    %first row is each Yposition (not binned, i think. 1986 in sample exp.)
    %each column after that is Fc3_DF for each neuron
    combinedYpos_no_rew=[ybinned(longruns_pos_no_rew) Fc3_DF(longruns_pos_no_rew,:)];
    
    combinedYneg_no_rew=[ybinned(longruns_neg_no_rew) Fc3_DF(longruns_neg_no_rew,:)];
    
     ginput(1);
    
    
    close all;
    
    
    
    
    
    

    %sorts ascending. I have no idea why this happens mwa
    sortedcombinedYpos_no_rew=sortrows(combinedYpos_no_rew);
   
    sortedcombinedYneg_no_rew=sortrows(combinedYneg_no_rew);
    
    %pos velocity only; no rewards
    pos_Fpos_no_rew=[];
    numbins=16; %orig 80
    sizebins=length(sortedcombinedYpos_no_rew)/numbins;
    numsigs_pos_no_rew=[];
    for j=1:numbins
        temp1=sortedcombinedYpos_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),1);
        temp1=temp1(find(~isnan(temp1)));
        pos_Fpos_no_rew(1,j)=mean(temp1);
        
        for jj=1:numneurons
            temp2=sortedcombinedYpos_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),jj+1);
            temp2=temp2(find(~isnan(temp2)));
            pos_Fpos_no_rew(jj+1,j)=mean(temp2);
            numsigs_pos_no_rew(jj,j)=sum(temp2>0)/length(temp2);
        end
    end
    
    
    %neg velocity only; no rewards
    pos_Fneg_no_rew=[];
    numbins=80;
    sizebins=length(sortedcombinedYneg_no_rew)/numbins;
    numsigs_neg_no_rew=[];
    for j=1:numbins     %j is each bin
        temp1=sortedcombinedYneg_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),1);
        temp1=temp1(find(~isnan(temp1)));
        pos_Fneg_no_rew(1,j)=mean(temp1);
        for jj=1:numneurons     %jj is each neuron
            temp2=sortedcombinedYneg_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),jj+1);
            temp2=temp2(find(~isnan(temp2)));
            pos_Fneg_no_rew(jj+1,j)=mean(temp2);
            numsigs_neg_no_rew(jj,j)=sum(temp2>0)/length(temp2);
        end
    end
    

    
    %%search for place fields in positive passes; no rewards
    tempypos_no_rew=pos_Fpos_no_rew(1,:);
    pos_rate_no_rew=[];
    fieldindpos_no_rew=zeros(numneurons,8);
    for j=1:numneurons
        fieldcount=1;
        temp1=smooth(pos_Fpos_no_rew(j+1,:),3);
        temp1sort=sort(temp1);
        base1=mean(temp1sort(1:round(baselength*numbins)));
        thresh2=(thresh1*(max(temp1)-mean(base1)))+base1;
        temp2=temp1>thresh2;
        L=bwlabeln(temp2,4);
        for ii=1:max(L)
            
            fieldchar(j).pos_no_rew(1,ii)=((max(temp1(find(L==ii)))));
            fieldchar(j).pos_no_rew(2,ii)=(tempypos_no_rew(max(find(L==ii)))-tempypos_no_rew(min(find(L==ii))));
            fieldchar(j).pos_no_rew(3,ii)=mean(numsigs_pos_no_rew(j,find(L==ii)));
            fieldchar(j).pos_no_rew(4,ii)=((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))));
            fieldchar(j).pos_no_rew(5,ii)=tempypos_no_rew(min(find(L==ii)));
            fieldchar(j).pos_no_rew(6,ii)=tempypos_no_rew(max(find(L==ii)));
            if ((max(temp1(find(L==ii)))))>minDF
                if (tempypos_no_rew(max(find(L==ii)))-tempypos_no_rew(min(find(L==ii))))>minfieldwidth
                    %if max(numsigs_pos_no_rew(j,find(L==ii)))>minrate
                    if mean(numsigs_pos_no_rew(j,find(L==ii)))>minrate
                        if ((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))))>minratio
                            fieldindpos_no_rew(j,fieldcount)=min(find(L==ii));
                            fieldindpos_no_rew(j,fieldcount+1)=max(find(L==ii))
                            fieldcount=fieldcount+2;
                            figure;
                            subplot(1,3,1)
                            hold on;
                            plot(tempypos_no_rew,temp1);
                            plot(tempypos_no_rew(min(find(L==ii)):max(find(L==ii))),temp1(min(find(L==ii)):max(find(L==ii))),'r');
                            title([num2str(j) ' pos no rew']);
                            xlim([-600 800]);
                            subplot(1,3,2)
                            stemtemp=combinedYpos_no_rew(:,j+1);
                            stemtemp1=find(stemtemp==0);
                            stemtemp(stemtemp1)=0.0075;
                            stem(combinedYpos_no_rew(:,1),stemtemp,'Markersize',.01);
                            title([num2str(j) ' pos no rew']);
                            xlim([-600 800]);
                            pos_rate_no_rew=[pos_rate_no_rew mean(numsigs_pos_no_rew(j,find(L==ii)))]
                            %suppress for large data sets.  121108  EBH
                            ginput(1);
                            %                             close all;
                        end
                    end
                end
            end
        end
    end
    
    
    
    
    %%search for place fields in negative passes; no rewards
    fieldindneg_no_rew=zeros(numneurons,8);
    tempyneg_no_rew=pos_Fneg_no_rew(1,:);
    neg_rate_no_rew=[];
    for j=1:numneurons
        fieldcount=1;
        temp1=smooth(pos_Fneg_no_rew(j+1,:),3);
        temp1sort=sort(temp1);
        base1=mean(temp1sort(1:round(.3*numbins)));
        thresh2=(thresh1*(max(temp1)-mean(base1)))+base1;
        temp2=temp1>thresh2;
        L=bwlabeln(temp2,4);
        for ii=1:max(L)
            
            fieldchar(j).neg_no_rew(1,ii)=((max(temp1(find(L==ii)))));
            fieldchar(j).neg_no_rew(2,ii)=(tempyneg_no_rew(max(find(L==ii)))-tempyneg_no_rew(min(find(L==ii))));
            fieldchar(j).neg_no_rew(3,ii)=mean(numsigs_neg_no_rew(j,find(L==ii)));
            fieldchar(j).neg_no_rew(4,ii)=((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))));
            fieldchar(j).neg_no_rew(5,ii)=tempyneg_no_rew(min(find(L==ii)));
            fieldchar(j).neg_no_rew(6,ii)=tempyneg_no_rew(max(find(L==ii)));
            if ((max(temp1(find(L==ii)))))>minDF
                if (tempyneg_no_rew(max(find(L==ii)))-tempyneg_no_rew(min(find(L==ii))))>minfieldwidth
                    %if max(numsigs_neg_no_rew(j,find(L==ii)))>minrate
                    if mean(numsigs_neg_no_rew(j,find(L==ii)))>minrate
                        if ((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))))>minratio
                            fieldindneg_no_rew(j,fieldcount)=min(find(L==ii));
                            fieldindneg_no_rew(j,fieldcount+1)=max(find(L==ii))
                            fieldcount=fieldcount+2;
                            figure; subplot(1,3,1)
                            hold on;
                            plot(tempyneg_no_rew,temp1);
                            plot(tempyneg_no_rew(min(find(L==ii)):max(find(L==ii))),temp1(min(find(L==ii)):max(find(L==ii))),'r');
                            title([num2str(j) ' neg no rew']);
                            xlim([-600 800]);
                            subplot(1,3,2);
                            stemtemp=combinedYneg_no_rew(:,j+1);
                            stemtemp1=find(stemtemp==0);
                            stemtemp(stemtemp1)=0.0075;
                            stem(combinedYneg_no_rew(:,1),stemtemp,'Markersize',.01);
                            title([num2str(j) ' neg no rew']);
                            xlim([-600 800]);
                            neg_rate_no_rew=[neg_rate_no_rew mean(numsigs_neg_no_rew(j,find(L==ii)))]
                            %suppress for large data sets.  121108  EBH
%                             ginput(1);
                            %                             close all;
                        end
                    end
                end
            end
        end
    end
    %
        ginput(1);
    %     ginput(1);
    %
    
%     close all
    
    
    %added below to get individual events and stats.  120921  EH
    IEI_str=('_IEI_');  %need to make string to concatenate for saving figs below
    Ave_str=('_Ave_');
    
    n_tot_ev=zeros(numneurons,1);%neuron total events
    Fc3_DF_events_label=zeros(length(Fc2),numneurons);    %preallocate
    
    for i=1:numneurons
        Fc3_DF_events_label(:,i)=bwlabel(Fc3_DF(:,i),4);
        % i
        n_tot_ev(i)=max(Fc3_DF_events_label(:,i));
        
        switch(n_tot_ev(i))
            case 0
                disp(i)
                disp(n_tot_ev(i))
                all_events(i)=struct;
                
            otherwise
                
                
                %all_events(i)=struct('info',zeros(n_tot_ev(i),12),'summary',zeros(2,12));
                all_events(i)=struct('info',zeros(12,(n_tot_ev(i))),'summary',zeros(12,2),'traces',zeros(236,n_tot_ev(i)),'ave',zeros(236,4));
                %236 frames
                all_events(i).info(11,n_tot_ev(i)) =NaN;    %sets last IEI to NaN
                
        end
    end
    
    
    for i= 1:numneurons
        switch(n_tot_ev(i))
            case 0
                
                
            otherwise
                
                for ev_loop=1:n_tot_ev(i)
                    all_events(i).info(1,ev_loop)=ev_loop;
                    all_events(i).info(2,ev_loop)=find(Fc3_DF_events_label(:,i)==ev_loop,1,'first');
                    all_events(i).info(3,ev_loop)=find(Fc3_DF_events_label(:,i)==ev_loop,1,'last');
                    all_events(i).info(4,ev_loop)= ((all_events(i).info(3,ev_loop))-(all_events(i).info(2,ev_loop)+1))*(ms_per_frame)/1000;
                    [m,p]=max(Fc3_DF(all_events(i).info(2,ev_loop):all_events(i).info(3,ev_loop),i));
                    all_events(i).info(5,ev_loop)=m;
                    all_events(i).info(6,ev_loop)=(p+(all_events(i).info(2,ev_loop))-1);
                    all_events(i).info(7,ev_loop)=ybinned(all_events(i).info(6,ev_loop));
                    %all_events(i).info(ev_loop,8)    %bin with max
                    all_events(i).info(9,ev_loop)=mean(Fc3_DF(all_events(i).info(2,ev_loop):all_events(i).info(3,ev_loop),i));
                    %          all_events(i).info(ev_loop,10)
                    %originally wanted row 12 for Instantaneous Freq for each event.
                    %Currently using it for overall freq (see below).
                    
                    
                    if (all_events(i).info(2,ev_loop))-20>1
                        if  (all_events(i).info(3,ev_loop)+215)<length(Fc2)
                            all_events(i).traces(:,ev_loop)=Fc2(((all_events(i).info(2,ev_loop))-20):((all_events(i).info(2,ev_loop))+215),i);
                        end
                    end
                    
                    
                end
                
                for ev_loop=1:(n_tot_ev(i)-1)
                    all_events(i).info(11,ev_loop)=(all_events(i).info(6,(ev_loop+1))-(all_events(i).info(6,ev_loop)))*(ms_per_frame/1000);
                    
                end
                all_events(i).traces(all_events(i).traces(:,:)==0)=NaN;
                for i_loop=2:12
                    all_events(i).summary(i_loop,1)=nanmean(all_events(i).info(i_loop,:));
                    all_events(i).summary(i_loop,2)=nanstd(all_events(i).info(i_loop,:));
                    
                end
                
                all_events(i).summary(1,1)=n_tot_ev(i,1);
                all_events(i).summary(12,1)=(n_tot_ev(i,1))/((length(Fc2)*ms_per_frame)/1000);
                %originally wanted row 12 for Instantaneous Freq for each event.
                %Currently using it for overall freq.
                
                for i_loop=1:236
                    all_events(i).ave(i_loop,1)=nanmean(all_events(i).traces(i_loop,:));
                    all_events(i).ave(i_loop,2)=nanstd(all_events(i).traces(i_loop,:));
                    
                end
                
                %          all_events(i).summary(2:12,:)=nanmean(all_events(i).info(2:12,:));
                %          all_events(i).summary(2:12,:)=nanstd(all_events(i).info(2:12,:));
                
                
                
                
                all_events(i).ave(:,3)=all_events(i).ave(:,1)+2*(all_events(i).ave(:,1));
                all_events(i).ave(:,4)=all_events(i).ave(:,1)-2*(all_events(i).ave(:,1));
                
                %figure;hold on; title(i); hist(all_events(i).info(11,:),0:((max(all_events(i).info(11,:)))/20):(max(all_events(i).info(11,:))));
                %figure;hold on; title(i); hist(all_events(i).info(11,:),0:0.5:80);
                
                figure;hold on; title(i); hist(all_events(i).info(11,:),0:0.25:6); xlim ([0 6]); %expand 0-4s int
                %To save IEI hist, uncomment line below
                %saveas (gcf,strcat((strcat(name, IEI_str)),int2str(i)),'emf');
                
                %figure;hold on; title(i); plot((all_events(i).traces(:,:)),'Color',[.8 .8 .8]);plot((all_events(i).ave(:,3)),':k');plot((all_events(i).ave(:,4)),':k');plot((all_events(i).ave(:,1)),'k','LineWidth',2);
                
                figure;hold on; title(i); plot(ms_per_frame:ms_per_frame:236*ms_per_frame,(all_events(i).traces(:,:)),'Color',[.8 .8 .8]);plot(ms_per_frame:ms_per_frame:236*ms_per_frame,(all_events(i).ave(:,3)),':k');plot(ms_per_frame:ms_per_frame:236*ms_per_frame,(all_events(i).ave(:,4)),':k');plot(ms_per_frame:ms_per_frame:236*ms_per_frame,(all_events(i).ave(:,1)),'k','LineWidth',2);
                %To save average transient, uncomment line below
                %saveas (gcf,strcat((strcat(name, Ave_str)),int2str(i)),'emf');
                %figure;plot(all_events(i).info(:,11));
                %figure;errorbar((all_events(i).ave(1,:))',(all_events(i).ave(2,:))');
                %suppress below for large data sets.  121108  EBH
%                 ginput(1);
                close all;
        end
    end
    
    pause(1);
    
    
    fieldsizes=[];  % this part just gets positive place field sizes and appends
    %to field sizes.  Doesn't keep track of which neuron it came from
    %or direction.
    
    for i=1:numneurons
%         if fieldindpos_rew(i,2)>1
%             %since it only looks at column 2 (place field#1 end bin), it
%             %will not include any cells with multiple place fields.
%             %If two fields pass above ifs, 2nd field will be in column
%             %3(start) and 4 (end) of fieldind.
%             %However, it's unlikely that this analysis will find any sig.
%             %double place fields since looking at the 25% lowest dF/F to
%             %define baseline will make it artificially high.
%             fieldsizes=[fieldsizes tempypos_rew(fieldindpos_rew(i,2))-tempypos_rew(fieldindpos_rew(i,1))];
%         end
        if fieldindpos_no_rew(i,2)>1
            fieldsizes=[fieldsizes tempypos_no_rew(fieldindpos_no_rew(i,2))-tempypos_no_rew(fieldindpos_no_rew(i,1))];
        end
%         if fieldindneg_rew(i,2)>1
%             fieldsizes=[fieldsizes tempyneg_rew(fieldindneg_rew(i,2))-tempyneg_rew(fieldindneg_rew(i,1))];
%         end
        if fieldindneg_no_rew(i,2)>1
            fieldsizes=[fieldsizes tempyneg_no_rew(fieldindneg_no_rew(i,2))-tempyneg_no_rew(fieldindneg_no_rew(i,1))];
        end
    end
    
    
    
    %convert to cm
    fieldsizes=fieldsizes*180/tracksize;
    % mean(fieldsizes)
    
    
    
    pos_placecells_no_rew=find(sum(fieldindpos_no_rew,2)>0);
    neg_placecells_no_rew=find(sum(fieldindneg_no_rew,2)>0);
    
    num_fields=sum(((sum(fieldindneg_no_rew,2)>0)+(sum(fieldindpos_no_rew,2)>0))>0);
    
    
    
    %    save(Fname,'fieldindpos_rew','fieldindpos_no_rew','fieldindneg_rew','fieldindneg_no_rew','pos_Fpos_no_rew','pos_Fpos_rew','pos_Fneg_no_rew','pos_Fneg_rew','minfieldwidth','minratio','thresh1','minDF','minrate','baselength','vthresh','ysize_thresh','trackstart','trackend','-append');
    %replaced Fname with Ffile below, untested.  120912  EH
    save(Ffile,'fieldindpos_rew','fieldindpos_no_rew','fieldindneg_rew','fieldindneg_no_rew','pos_Fpos_no_rew','pos_Fpos_rew','pos_Fneg_no_rew','pos_Fneg_rew','minfieldwidth','minratio','thresh1','minDF','minrate','baselength','vthresh','ysize_thresh','trackstart','trackend','-append');
    
    
    
    
    if num_fields>0
        
        
        %shuffle test
        
        %%do shuffle test
        numrand=1000;
%         total_pos_randfields_rew=zeros(numneurons,1);
        total_pos_randfields_no_rew=zeros(numneurons,1);
%         total_neg_randfields_rew=zeros(numneurons,1);
        total_neg_randfields_no_rew=zeros(numneurons,1);
        
        for hh=1:numrand
            hh
%             combinedYpos_rew=[];
            combinedYpos_no_rew=[];
%             combinedYneg_rew=[];
            combinedYneg_no_rew=[];
            
            %randomly shuffle Fc3
            new_Fc3_DF=[];
            for ii=1:numneurons
                
                %find starts and ends of transients in Fc3 files
                ends=find(diff(Fc3_DF(:,ii)>0)<0)+1;    %end is actually one frame after end
                starts=find(diff(Fc3_DF(:,ii)>0)>0);    %start is one frame before start
                
                %deal with end cases
                if size(starts,1)>size(ends,1)
                    ends=[ends' length(Fc3_DF(:,ii))];
                end
                if size(starts,1)<size(ends,1)
                    starts=[0 starts']; % why ' ????
                    %starts=[1 starts'];    %this is used for subdivided F
                    %files.  Not sure why 0 doesn't work.  See notes in
                    %_subdivided shuffle for more detail
                    
                end
                
                trans_dur=0;
                for i=1:size(starts,1)
                    trans_dur=trans_dur+(ends(i)-starts(i));
                end
                mean_trans_dur=round(trans_dur/size(starts,1));
                num_randbins=round(length(Fc3_DF(:,ii))/mean_trans_dur);
                size_randbins=round(length(Fc3_DF(:,ii))/num_randbins);
                
                %to insure num_rand<<possibilities (num_rand=10000)
                if or(isnan(num_randbins),num_randbins<9)
                    num_randbins=9;
                    size_randbins=round(length(Fc3_DF(:,ii))/num_randbins);
                end
                
                bin_rand=randperm(num_randbins);
                new_Fc3=0;
                for j=1:num_randbins
                    if bin_rand(j)<=size(starts,1)
                        new_Fc3=[new_Fc3 Fc3_DF(starts(bin_rand(j)):ends(bin_rand(j)),ii)'];
                    else
                        new_Fc3=[new_Fc3 zeros(size_randbins,1)'];
                    end
                    
                end
                
                if size(new_Fc3,2)>length(Fc3_DF(:,ii))
                    new_Fc3=new_Fc3(2:length(Fc3_DF(:,ii))+1);
                end
                %if size of new_Fc3 is too small, insert zeros at random location
                if size(new_Fc3,2)<length(Fc3_DF(:,ii))
                    zer_pos=round(rand*length(new_Fc3));
                    size_diff=length(Fc3_DF(:,ii))-size(new_Fc3,2);
                    new_Fc3=[new_Fc3(1:zer_pos) zeros(size_diff,1)' new_Fc3(zer_pos+1:end)];
                end
                new_Fc3_DF=[new_Fc3_DF;new_Fc3];
            end
            
            new_Fc3_DF=new_Fc3_DF';
            %bin position and average Ftrace
            combinedYpos_rew=[ybinned(longruns_pos_rew) new_Fc3_DF(longruns_pos_rew,:)];
            combinedYpos_no_rew=[ybinned(longruns_pos_no_rew) new_Fc3_DF(longruns_pos_no_rew,:)];
            combinedYneg_rew=[ybinned(longruns_neg_rew) new_Fc3_DF(longruns_neg_rew,:)];
            combinedYneg_no_rew=[ybinned(longruns_neg_no_rew) new_Fc3_DF(longruns_neg_no_rew,:)];
            sortedcombinedYpos_rew=sortrows(combinedYpos_rew);
            sortedcombinedYpos_no_rew=sortrows(combinedYpos_no_rew);
            sortedcombinedYneg_rew=sortrows(combinedYneg_rew);
            sortedcombinedYneg_no_rew=sortrows(combinedYneg_no_rew);
            
            
            %pos velocity only; rewards
            pos_Fpos_rew=[];
            numbins=80;
            sizebins=floor(length(sortedcombinedYpos_rew)/numbins);
            numsigs_pos_rew=[];
            for j=1:numbins
                temp1=sortedcombinedYpos_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),1);
                temp1=temp1(find(~isnan(temp1)));
                pos_Fpos_rew(1,j)=mean(temp1);
                for jj=1:numneurons
                    temp2=sortedcombinedYpos_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),jj+1);
                    temp2=temp2(find(~isnan(temp2)));
                    pos_Fpos_rew(jj+1,j)=mean(temp2);
                    numsigs_pos_rew(jj,j)=sum(temp2>0)/length(temp2);
                end
            end
            
            %pos velocity only; no rewards
            pos_Fpos_no_rew=[];
            numbins=80;
            sizebins=floor(length(sortedcombinedYpos_no_rew)/numbins);
            numsigs_pos_no_rew=[];
            for j=1:numbins
                temp1=sortedcombinedYpos_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),1);
                temp1=temp1(find(~isnan(temp1)));
                pos_Fpos_no_rew(1,j)=mean(temp1);
                for jj=1:numneurons
                    temp2=sortedcombinedYpos_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),jj+1);
                    temp2=temp2(find(~isnan(temp2)));
                    pos_Fpos_no_rew(jj+1,j)=mean(temp2);
                    numsigs_pos_no_rew(jj,j)=sum(temp2>0)/length(temp2);
                end
            end
            
            
            
            %neg velocity only; rewards
            pos_Fneg_rew=[];
            numbins=80;
            sizebins=floor(length(sortedcombinedYneg_rew)/numbins);
            numsigs_neg_rew=[];
            for j=1:numbins
                temp1=sortedcombinedYneg_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),1);
                temp1=temp1(find(~isnan(temp1)));
                pos_Fneg_rew(1,j)=mean(temp1);
                for jj=1:numneurons
                    temp2=sortedcombinedYneg_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),jj+1);
                    temp2=temp2(find(~isnan(temp2)));
                    pos_Fneg_rew(jj+1,j)=mean(temp2);
                    numsigs_neg_rew(jj,j)=sum(temp2>0)/length(temp2);
                end
            end
            
            
            %neg velocity only; no rewards
            pos_Fneg_no_rew=[];
            numbins=80;
            sizebins=floor(length(sortedcombinedYneg_no_rew)/numbins);
            numsigs_neg_no_rew=[];
            for j=1:numbins
                temp1=sortedcombinedYneg_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),1);
                temp1=temp1(find(~isnan(temp1)));
                pos_Fneg_no_rew(1,j)=mean(temp1);
                for jj=1:numneurons
                    temp2=sortedcombinedYneg_no_rew(round((j-1)*sizebins+1):round(floor(j*sizebins)),jj+1);
                    temp2=temp2(find(~isnan(temp2)));
                    pos_Fneg_no_rew(jj+1,j)=mean(temp2);
                    numsigs_neg_no_rew(jj,j)=sum(temp2>0)/length(temp2);
                end
            end
            
            
            
            %%search for place fields in positive passes; rewards
            fieldindpos_rew=zeros(numneurons,8);
            tempypos_rew=pos_Fpos_rew(1,:);
            for j=1:numneurons
                fieldcount=1;
                temp1=smooth(pos_Fpos_rew(j+1,:),3);
                temp1sort=sort(temp1);
                base1=mean(temp1sort(1:round(baselength*numbins)));
                thresh2=(thresh1*(max(temp1)-mean(base1)))+base1;
                temp2=temp1>thresh2;
                L=bwlabeln(temp2,4);
                for ii=1:max(L)
                    if ((max(temp1(find(L==ii)))))>minDF
                        if (tempypos_rew(max(find(L==ii)))-tempypos_rew(min(find(L==ii))))>minfieldwidth
                            %if max(numsigs_pos_rew(j,find(L==ii)))>minrate
                            if mean(numsigs_pos_rew(j,find(L==ii)))>minrate
                                if ((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))))>minratio
                                    fieldindpos_rew(j,fieldcount)=min(find(L==ii));
                                    fieldindpos_rew(j,fieldcount+1)=max(find(L==ii));
                                    fieldcount=fieldcount+2;
                                    total_pos_randfields_rew(j)=total_pos_randfields_rew(j)+1;
                                end
                            end
                        end
                    end
                end
            end
            
            %%search for place fields in positive passes; no rewards
            fieldindpos_no_rew=zeros(numneurons,8);
            tempypos_no_rew=pos_Fpos_no_rew(1,:);
            for j=1:numneurons
                fieldcount=1;
                temp1=smooth(pos_Fpos_no_rew(j+1,:),3);
                temp1sort=sort(temp1);
                base1=mean(temp1sort(1:round(baselength*numbins)));
                thresh2=(thresh1*(max(temp1)-mean(base1)))+base1;
                temp2=temp1>thresh2;
                L=bwlabeln(temp2,4);
                for ii=1:max(L)
                    if ((max(temp1(find(L==ii)))))>minDF
                        if (tempypos_no_rew(max(find(L==ii)))-tempypos_no_rew(min(find(L==ii))))>minfieldwidth
                            %if max(numsigs_pos_no_rew(j,find(L==ii)))>minrate
                            if mean(numsigs_pos_no_rew(j,find(L==ii)))>minrate
                                if ((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))))>minratio
                                    fieldindpos_no_rew(j,fieldcount)=min(find(L==ii));
                                    fieldindpos_no_rew(j,fieldcount+1)=max(find(L==ii));
                                    fieldcount=fieldcount+2;
                                    total_pos_randfields_no_rew(j)=total_pos_randfields_no_rew(j)+1;
                                end
                            end
                        end
                    end
                end
            end
            
            
            %%search for place fields in negative passes; rewards
            fieldindneg_rew=zeros(numneurons,8);
            tempyneg_rew=pos_Fneg_rew(1,:);
            for j=1:numneurons
                fieldcount=1;
                temp1=smooth(pos_Fneg_rew(j+1,:),3);
                temp1sort=sort(temp1);
                base1=mean(temp1sort(1:round(baselength*numbins)));
                thresh2=(thresh1*(max(temp1)-mean(base1)))+base1;
                temp2=temp1>thresh2;
                L=bwlabeln(temp2,4);
                for ii=1:max(L)
                    if ((max(temp1(find(L==ii)))))>minDF
                        if (tempyneg_rew(max(find(L==ii)))-tempyneg_rew(min(find(L==ii))))>minfieldwidth
                            %if max(numsigs_neg_rew(j,find(L==ii)))>minrate
                            if mean(numsigs_neg_rew(j,find(L==ii)))>minrate
                                if ((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))))>minratio
                                    fieldindneg_rew(j,fieldcount)=min(find(L==ii));
                                    fieldindneg_rew(j,fieldcount+1)=max(find(L==ii));
                                    fieldcount=fieldcount+2;
                                    total_neg_randfields_rew(j)=total_neg_randfields_rew(j)+1;
                                end
                            end
                        end
                    end
                end
            end
            
            
            %%search for place fields in negative passes; no
            %%rewards
            fieldindneg_no_rew=zeros(numneurons,8);
            tempyneg_no_rew=pos_Fneg_no_rew(1,:);
            for j=1:numneurons
                fieldcount=1;
                temp1=smooth(pos_Fneg_no_rew(j+1,:),3);
                temp1sort=sort(temp1);
                base1=mean(temp1sort(1:round(baselength*numbins)));
                thresh2=(thresh1*(max(temp1)-mean(base1)))+base1;
                temp2=temp1>thresh2;
                L=bwlabeln(temp2,4);
                for ii=1:max(L)
                    if ((max(temp1(find(L==ii)))))>minDF
                        if (tempyneg_no_rew(max(find(L==ii)))-tempyneg_no_rew(min(find(L==ii))))>minfieldwidth
                            %if max(numsigs_neg_no_rew(j,find(L==ii)))>minrate
                            if mean(numsigs_neg_no_rew(j,find(L==ii)))>minrate
                                if ((mean(temp1(find(L==ii))))/mean(temp1(find(L~=ii))))>minratio
                                    fieldindneg_no_rew(j,fieldcount)=min(find(L==ii));
                                    fieldindneg_no_rew(j,fieldcount+1)=max(find(L==ii));
                                    fieldcount=fieldcount+2;
                                    total_neg_randfields_no_rew(j)=total_neg_randfields_no_rew(j)+1;
                                end
                            end
                        end
                    end
                end
                
                
            end
            
            
        end
        
        
        
        aa_rew=total_pos_randfields_rew/numrand;
        aa_no_rew=total_pos_randfields_no_rew/numrand;
        bb_rew=total_neg_randfields_rew/numrand;
        bb_no_rew=total_neg_randfields_no_rew/numrand;
        positives_rew=[pos_placecells_rew aa_rew(pos_placecells_rew)]
        positives_no_rew=[pos_placecells_no_rew aa_no_rew(pos_placecells_no_rew)]
        negatives_rew=[neg_placecells_rew bb_rew(neg_placecells_rew)]
        negatives_no_rew=[neg_placecells_no_rew bb_no_rew(neg_placecells_no_rew)]
        
        num_repeats=0;
        for_repeats=[pos_placecells_rew;pos_placecells_no_rew;neg_placecells_rew;neg_placecells_no_rew];
        
        for tt=1:max(for_repeats)
            if length(find(for_repeats==tt))>1
                num_repeats=num_repeats+1;
            end
        end
        maxp=0.05;
        numplacecells=length(pos_placecells_rew)+length(pos_placecells_no_rew)+length(neg_placecells_rew)+length(neg_placecells_no_rew)-num_repeats;
        sig_place_cells=([aa_rew(pos_placecells_rew);aa_no_rew(pos_placecells_no_rew);bb_rew(neg_placecells_rew);bb_no_rew(neg_placecells_no_rew)])<maxp;
        rates=[pos_rate_rew pos_rate_no_rew neg_rate_rew neg_rate_no_rew];
        sig_rates=rates(sig_place_cells);
        sig_sizes=fieldsizes(sig_place_cells);
        listcells=[pos_placecells_rew;pos_placecells_no_rew;neg_placecells_rew;neg_placecells_no_rew];
        sig_placecell_numbers=listcells(sig_place_cells);
        sig_placecell_numbers_positive=0;
        sig_place_cells_positive=([aa_rew(pos_placecells_rew);aa_no_rew(pos_placecells_no_rew)])<maxp;
        list_positive=[pos_placecells_rew;pos_placecells_no_rew];
        sig_placecell_numbers_positive=list_positive(sig_place_cells_positive);
        sig_placecell_numbers_neg=0;
        sig_place_cells_neg=([bb_rew(neg_placecells_rew);bb_no_rew(neg_placecells_no_rew)])<maxp;
        list_neg=[neg_placecells_rew;neg_placecells_no_rew];
        sig_placecell_numbers_neg=list_neg(sig_place_cells_neg);
        
        
        %look if place cell has multiple fields (in diff dir or no rew, not
        %multipeaked).
        num_sig_repeats=0;
        for tt=1:max(sig_placecell_numbers)
            if length(find(sig_placecell_numbers==tt))>1
                num_sig_repeats=num_sig_repeats+1;
            end
        end
        
        num_sig_placecells=length(sig_placecell_numbers)-num_sig_repeats
        %only counted as place cell once
    else
        sig_placecell_numbers=[];
        sig_placecell_numbers_neg=[];
        sig_placecell_numbers_positive=[];
        
        maxp=0.05;
        numrand=[];
        positives_rew=[];
        positives_no_rew=[];
        negatives_rew=[];
        negatives_no_rew=[];
        numplacecells=[];
        aa_rew=[];
        aa_no_rew=[];
        bb_rew=[];
        bb_no_rew=[];
        sig_rates=[];
        sig_sizes=[];
        rates=[];
        num_sig_placecells=[];
    end
    
    
    
    % save(Fname,'num_sig_placecells','sig_placecell_numbers','maxp','numrand','positives_rew','positives_no_rew','negatives_rew','negatives_no_rew','numplacecells','aa_rew','aa_no_rew','bb_rew','bb_no_rew','pos_rate_rew','pos_rate_no_rew','neg_rate_rew','neg_rate_no_rew','sig_rates','sig_sizes','rates','fieldsizes','-append');
    %replaced Fname with Ffile below.   120912  EH
    save(Ffile,'Fc3_DF_events_label','n_tot_ev','all_events','fieldindpos_rew','fieldindpos_no_rew','fieldindneg_rew','fieldindneg_no_rew','pos_Fpos_no_rew','pos_Fpos_rew','pos_Fneg_no_rew','pos_Fneg_rew',...
        'minfieldwidth','minratio','thresh1','minDF','minrate','baselength','vthresh','ysize_thresh','trackstart','trackend','num_sig_placecells','sig_placecell_numbers','pos_rate_rew','pos_rate_no_rew',...
        'neg_rate_rew','neg_rate_no_rew','fieldsizes','fieldchar','sig_placecell_numbers_neg','sig_placecell_numbers_positive','-append');
    
    
    
    
    
    
end









