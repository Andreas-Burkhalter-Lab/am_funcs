function[type]=what_type_one_line(wis,programming)
%determine the difference between complex, upwards, downwards,flat and
%chervron
%1.Complex calls displayed one syllable containing two or more directional changes in pitch, each ?6.25 kHz.
%DONE%4.Upward-modulated calls exhibited a continuous increase in pitch that was ?12.5 kHz, with a terminal dominant frequency at least 6.25 kHz more than the pitch at the beginning of the vocalization.
%DONE%5.Downward-modulated calls exhibited a continuous decrease in pitch that was ?12.5 kHz, with a terminal dominant frequency at least 6.25 kHz less than the pitch at the beginning of the vocalization.
%DONE 10.Flat calls displayed a constant beginning and the ending of the pitch frequency remained constant (?3 kHz of each other).
%6.Chevron calls resembled an ‘inverted-U’, which was identified by a continuous increase in pitch ?12.5 kHz followed by a decrease that was ?6.25 kHz.


[a pit] = max(wis);

pit(pit==1 | pit==0)=[];%bug proof  ... AM combined into one OR statement


%first take out extremes remember at this point they are supposed to be
%basically sine waves %%this is 25 hrtz in one direction in 1 millescond
%maybe it should be smaller
%figure;plot(pit,'b*');
% for x=1:length(pit)-1
%     if (abs(pit(x)-pit(x+1))>50) ~=1
%         newpit(x)=pit(x);
%     end
% end
% pit=[newpit pit(end)];%awkward;
% pit(pit==0)=[];


%smooth it
%hold on;plot(pit,'r');
if size(pit,2)>4
spit=smooth(pit,5);
spit(1)=pit(1);
spit(2) = mean(pit(1:2));
spit(3) = mean(pit(1:3));
spit(4) = mean(pit(1:4));
else 
spit=pit;
end
%hold on;plot(spit,'b');
%direction=diff(spit);

%get local maximinums and minimums
%hving trouble with little fluxes, let smooth again
dif=diff(spit);
if length(spit)<5
    type=17;
    return
end
sdif=smooth(dif,3);
sdif(2)=mean(dif(1:2));
xmax = find(diff(sdif<.001)==1)+1;
xmin = find(diff(sdif>.001)==1)+1;%this helps with plateus not sure why I need it
idx=ones(size(spit));

size_consec_change=diff(spit(sort([1 xmax xmin length(spit)])));

if 1
    idx(1)=143;
    idx(end)=143;
    idx(xmin)=repmat(140,size(xmin));
    idx(xmax)=repmat(145,size(xmax));
    idx(xmin)=repmat(140,size(xmin));
    figure;plot(spit);hold on;plot(idx,'k*');
end

%%flat%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flat calls displayed a constant beginning and the ending of the pitch frequency
%remained constant (3 kHz of each other(6 in half freqs)).
therange = max(pit)-min(pit);
if therange<=3*2%half freq steps
    type=10;
    return
end

%%%%%%Complex%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Complex calls displayed one syllable containing two
%or more directional changes in pitch, each ?6.25 kHz.
pos= length(find(size_consec_change>6.25*2));
neg= length(find(size_consec_change<-6.25*2));
if pos && neg
    %%chervron
    %7.Chevron calls resembled an ‘inverted-U’, which was identified by a continuous increase in pitch ?12.5 kHz followed by a decrease that was ?6.25 kHz.
    if sum(size_consec_change>=12.5*2)%12.5khx*2 for half steps
        idf = find(size_consec_change>=12.5*2);
        if sum(size_consec_change(idf:end)<=-6.25*2)==1
            type=6;
        end
        type=1;
    end
end
%%%%%%%%%%%%%%almost chevron
pos= length(find(size_consec_change>10));
neg= length(find(size_consec_change<-10));
if pos && neg
    %%chervron
    %7.Chevron calls resembled an ‘inverted-U’, which was identified by a continuous increase in pitch ?12.5 kHz followed by a decrease that was ?6.25 kHz.
    if sum(size_consec_change>=10)%12.5khx*2 for half steps
        idf = find(size_consec_change>=10);
        if sum(size_consec_change(idf:end)<=-10)==1
            type=16;
            return
        end
    end
end


%Downwards%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Downward-modulated calls exhibited a continuous decrease in pitch that was ?12.5 kHz, with a terminal
%dominant frequency at least 6.25 kHz less than the pitch at the beginning of the vocalization.
if pit(1)-pit(end)>6.25*2%with a terminal%dominant frequency at least 6.25 kHz less than the pitch at the beginning of the vocalization.
    if sum(size_consec_change<=-12.5*2)%12.5khx*2 for half steps
        type=5;
        return
    end
end

%%%Upwards%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Upward-modulated calls exhibited a continuous increase in pitch that was ?12.5 kHz, with a terminal
%  dominant frequency at least 6.25 kHz more than the pitch at the beginning of the vocalization.
if  pit(end)-pit(1)>6.25*2
    if sum(size_consec_change>=12.5*2)%12.5khx*2 for half steps
        type=4;
        return
    end
end


%%%%almost down
if sum(size_consec_change<=-6.25*2)
    type=15;
    return
end

%almost up
%%%%almost down
if sum(size_consec_change>=6.25*2)
    type=14;
    return
end

%%%%%almost flat
if therange<=20
    type=20;
    return
end

figure;subplot(2,2,1);imagesc(wis);hold on;  set(gca,'YDir','normal');set(gca,'yticklabel',(0:25:250));
subplot(2,2,2);plot(spit);hold on; plot(idx,'*r');
type=99;
disp(type);