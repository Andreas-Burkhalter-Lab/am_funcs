function [  syl_type]= scattoni(whis)
%from scattoni et al.  <whis> is n by 1 cell array.  Each cell has the
%whistle (whis from tims program)
%1.	Complex calls displayed one sound component containing two or more directional changes in frequency with varying direction and rate of changes, each ?6.25 kHz.
%2.	Unstructured calls did not contain a main sound component and their shape cannot be assimilated to any of the other categories (see figure S1 for details).
%3.	Two-component calls consisted of two components: a main call (with a flat or downward frequency change) with an additional short component of higher frequency.
%4.	Upward-modulated calls exhibited a continuous increase in frequency that was ?12.5 kHz, with a terminal dominant frequency at least 6.25 kHz more than the frequency at the beginning of the vocalization.
%5.	Downward-modulated calls exhibited a continuous decrease in frequency that was ?12.5 kHz, with a terminal dominant frequency at least 6.25 kHz less than the initial frequency of the vocalization.
%6.	Chevron calls resembled an ‘inverted-U shape’, which was identified by a continuous increase in frequency ?12.5 kHz followed by a decrease that was ?6.25 kHz.
%7.	Short calls were punctuated and shorter than 5 ms.
%8.	Composite calls/harmonics were formed by two parallel lines consisting of a continuous sound with a minor frequency modulation and a single higher amplified harmonic component.
%9.	Frequency steps calls consisted of three components: a main call resembling the composite calls flanked by two discontinuous “step” on the sonographic display, but with no gaps on the time scale.
%10.	Flat calls presented a constant frequency, including the initial and terminal frequency with changes ?3 kHz.
%11.	Other 2< lines
%12.
%13.
%14.	14 almost up
%15.	15 Almost down
%16.	 A
%17.	A
%18.	Almost composite (overlap 50 percent)
%19.	L
%20.	Almost flat
%set to zero;

programming=0;
judge=0;

max_short_duration = 5;%in millseconds
sizeofnoise = 2;%how many neighbors have to have a signal for a pixel not to be considered noise
%if default (1) only isolated pixels will be ignored.
per_of_max = .08; %threshold of noise




for syl = 1:size(whis,2)
    close all
    wis=full(abs(whis{syl}));
    
    %SHORT SYLLABLE
    [a pit] = max(wis);
    pit(pit==1 | pit==0)=[];%bug proof  ... AM combined into one OR statement
    
    if size(pit,2)< max_short_duration+1
        syl_type(syl)=7;
    else%ELSE
        % 1 line or more?
        [wis numlines(syl) CC]=determine_line_number(syl,wis,sizeofnoise,per_of_max,programming);
        
        %more than one line? %automatic classfication fork
        if numlines(syl)==1
            syl_type(syl) = what_type_one_line(wis,judge);
            %  imagesc(wis);   set(gca,'YDir','normal');set(gca,'yticklabel',(0:25:250));
        else
            syl_type(syl) = what_type_morethanone_one(wis,numlines(syl),CC);
            
        end
        
        
    end
end
if programming
    rightanswers=xlsread('d:/numlines.xlsx');
    grader=abs(numlines'-rightanswers(1:size(numlines,2),1));
    grader(grader==.5)=0;
    sprintf('your grade is %d',sum(grader))
    sprintf('with overestimation in %d',length(find(numlines'-rightanswers(1:size(numlines,2),1)>.6)))
    sprintf('with underestimation in %d',length(find(numlines'-rightanswers(1:size(numlines,2),1)<-.6)))
    for x=1:size(grader,1)
        if grader(x)<.6
            close(x);
        end
    end
end
end



function  [wis_g numlines CC] = determine_line_number(syl,wis,sizeofnoise,per_of_max,programming)
if programming
    figure(syl)
    hold on;
    subplot(2,8,1)
    imagesc(wis)
    set(gca,'YDir','normal');
    colormap(flipud(gray));
    
    
    
    set(gca,'YTickLabel',{'25','50','75','100','125'});
    title(syl);
    set(gca,'XTick',size(wis,2));
    if rem(syl,24)==0
        figure;
        hold on;
        su = 1;
    end
end
%get rid of any single pixels
%first temparary added zeros to each side
wis_t=zeros(size(wis,1)+2,size(wis,2)+2);
wis_t(2:end-1,2:end-1)=wis;
%i switched row and col
for c=2:size(wis_t,2)-1
    for row=2:size(wis_t,1)-1
        nn = [ wis_t(row+1,c),wis_t(row-1,c),wis_t(row+1,c+1),wis_t(row+1,c-1),wis_t(row,c+1),...
            wis_t(row,c-1),wis_t(row-1,c-1),wis_t(row-1,c+1)];
        
        if   length(find(nn))< sizeofnoise
            wis_t(row,c)=0;
        end
    end
end
wis=wis_t(2:end-1,2:end-1);
wis_g=wis;
%next let's smooth
b=ones(1,5)/5;
a=1;
clear filteredData;
for col=1:size(wis,2)
    filteredData(:,col)=filter(b,a,wis(:,col));
end
wis=filteredData;
clear filteredData;
for row=1:size(wis,1)
    filteredData(row,:)=filter(b,a,wis(row,:));
end
wis=filteredData;
clear filteredData;

if programming
    subplot(2,8,2)
    imagesc(wis)
    set(gca,'YDir','normal');
    colormap(flipud(gray));
    %first lets filter
end
[themax max_loc] =(max(wis));
clear highs;
for x=1:size(wis,2)
    highs(size(wis,1),x)=zeros;
    highs(find(wis(:,x)>per_of_max*mean(themax)),x)=1;
end
if programming
    subplot(2,8,3)
    imagesc(highs);
    set(gca,'YDir','normal');
    colormap(flipud(gray));
    wis=highs;
end

%smooth again

for col=1:size(wis,2)
    filteredData(:,col)=filter(b,a,wis(:,col));
end
wis=filteredData;
clear filteredData;

b=ones(1,7)/7;
a=1;
for row=1:size(wis,1)
    filteredData(row,:)=filter(b,a,wis(row,:));
end
wis=filteredData;
clear filteredData;

if programming
    subplot(2,8,4)
    imagesc(wis);
    set(gca,'YDir','normal');
    colormap(flipud(gray));
    hold on;
end

CC = bwconncomp(wis);
clear freqlines;
for x=1:length(CC.PixelIdxList)
    if (size(CC.PixelIdxList{x},1))>4
        freqlines(x)=1;
    end
end
if isempty(CC.PixelIdxList)
    freqlines=0
else
    freqlines=sum(freqlines);
    if freqlines==0
        freqlines=1;
    end
end
numlines = freqlines;

if  programming
    subplot(2,8,4)
    title(sprintf('%2.0d',freqlines))
    
end


end

