function [posL,pos,negL,neg] = split_paths(ybinned,varargin)
%returns posL, negL - labeled blocks above run threshold
%pos, neg - indices of these runs
if nargin>1
    Fs=varargin{1};
    cropped_length=length(ybinned)-mod(length(ybinned),round(Fs/4));
    nybinned=interp1(1:cropped_length,smooth(ybinned(1:cropped_length)),...
        linspace(1,cropped_length,cropped_length/(Fs/4)))';
    nybinned=decimate(ybinned(1:cropped_length),round(Fs/4));
nybinned=smooth(nybinned);
else
    nybinned=ybinned;
end

ysize_thresh=peak2peak(ybinned)/4;%Have to go at least 40cm
vthresh=(ysize_thresh/20)/Fs; %need to be going at least 5cm/s
% figure;
% plot(ybinned);

ydiff=diff(nybinned);
if size(ydiff>0)

ydiff=[ydiff(1); ydiff];



%positive direction

L = bwlabel(ydiff>vthresh); %120508 EH
for i=1:max(L)
    ysize=abs(nybinned(find(L==i,1,'first'))-nybinned(find(L==i,1,'last')));
    if ysize<ysize_thresh
        L(L==i)=0;
    end
end
posL=L;
pos=find(L>0);


%negative direction

L = bwlabel(ydiff<-vthresh);    %120508 EH
for i=1:max(L)
    ysize=abs(nybinned(find(L==i,1,'last'))-nybinned(find(L==i,1,'first')));
    if ysize<ysize_thresh
        L(L==i)=0;
    end
end
negL=L;
neg=find(L>0);


if nargin>1
% posL=repelem(posL,round(length(ybinned)/length(nybinned)));
posL=interp1(1:length(nybinned),posL',linspace(1,length(nybinned),length(ybinned)));
% posL=[posL; zeros(mod(length(ybinned),round(length(ybinned)/length(nybinned))),1)];

% negL=repelem(negL,round(length(ybinned)/length(nybinned)));
negL=interp1(1:length(nybinned),negL,linspace(1,length(nybinned),length(ybinned)));
% negL=[negL; zeros(mod(length(ybinned),round(length(ybinned)/length(nybinned))),1)];
pos=find(posL>0);
neg=find(negL>0);


end
else 
    pos=nan; neg=nan; posL=nan; negL=nan;
 end   
end

