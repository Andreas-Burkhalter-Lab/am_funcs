
function Fc3_DF=calc_Fc3(Fc2)
find_min_signif_dur_Fpoly
%define baseline, find start, duration and max std amplitude of
%transients in both positive and neg. direction.
for i=1:numtraces
    baselineind=find(Fc2(:,i)<upperbase(i));
    baseline=Fc2(baselineind,i);
    basemedian=median(baseline);
    basestd=std(baseline);
    Fc2std(:,i)=(Fc2(:,i)-basemedian)./basestd;
    %positive going transients
    Fc2state(:,i)=double_thresh(Fc2std(:,i),2.0,0.5);
    upticks=find(diff(Fc2state(:,i))==1);
    downticks=find(diff(Fc2state(:,i))==-1);
    %deal with end cases
    if numel(upticks)>0
        if upticks(1)>downticks(1)
            upticks=[1 upticks'];
        end
        
        if upticks(end)>downticks(end)
            downticks=[downticks' length(Fc2)];
        end
    else
        transdur(i)=0;
        transmax(i)=0;
        
    end
    if length(upticks)>0
        for j=1:length(upticks)
            transmax(i,j)=max(Fc2std(upticks(j):downticks(j),i));
            transdur(i,j)=downticks(j)-upticks(j);
            transupticks(i,j)=upticks(j);
        end
    else
        transmax(i,j)=0;
        transdur(i,j)=0;
        transupticks(i,j)=0;
    end
    
    
    %negative going transients
    negFc2state(:,i)=double_thresh(-Fc2std(:,i),2.0,0.5);
    negupticks=find(diff(negFc2state(:,i))==1);
    negdownticks=find(diff(negFc2state(:,i))==-1);
    if numel(negupticks)>0
        if negupticks(1)>negdownticks(1)
            negupticks=[1 negupticks'];
        end
        if negupticks(end)>negdownticks(end)
            negdownticks=[negdownticks' length(Fc2)];
        end
    else
        negtransmax(i)=0;
        negtransmax(i)=0;
    end
    if length(negupticks)>0
        for j=1:length(negupticks)
            negtransmax(i,j)=max(-Fc2std(negupticks(j):negdownticks(j),i));
            negtransdur(i,j)=negdownticks(j)-negupticks(j);
            negtransupticks(i,j)=negupticks(j);
        end
    else
        negtransmax(i,j)=0;
        negtransdur(i,j)=0;
        negtransupticks(i,j)=0;
    end
    
end




%Make Fc3 files by only plotting transients <5% error rate and all else
%0; and only plotting std>3
% figure;
for i=1:numtraces
%     disp(['Currently on Cell: ',num2str(i)]);
    std_sig=zeros(length(Fc2),5);
%     clf;
%     hold on;
    
    for k=3:5
        
        nnn=0;
        if k==5
            nnn=(transmax(i,:)>k).*transdur(i,:);
            nnn=nnn.*(nnn>std_sig_dur_min_day(k));
        else
            nnn=(and(transmax(i,:)<(k+1),transmax(i,:)>k)).*transdur(i,:);
            nnn=nnn.*(nnn>std_sig_dur_min_day(k));
        end
        
        for jj=1:length(transupticks(i,:))
            if(transupticks(i,jj)>0)
                if(nnn(jj)>0)
                    std_sig((transupticks(i,jj)+1):(transupticks(i,jj)+nnn(jj)),k)=1;
                end
            end
        end
    end
    
    Fc3(:,i)=sum(std_sig,2).*Fc2std(:,i);
%     plot(Fc2std(:,i));
%     hold on;
%     plot(sum(std_sig,2).*Fc2std(:,i),'r');
%     ginput(1);
%     hold off;
    
end
% ginput(1);
%Save Fc3 in F file


%turn Fc3 back into DF/F (from std)
Fc3_DF=zeros(size(Fc2));
for i=1:numtraces
    baselineind=find(Fc2(:,i)<upperbase(i));
    baseline=Fc2(baselineind,i);
    basemedian=median(baseline);
    basestd=std(baseline);
    Fc3_DF(:,i)=Fc3(:,i)*basestd;
    
end

% save(fullFfile,'Fc3','Fc3_DF','-append');
% 
% clear all;
% close all;

