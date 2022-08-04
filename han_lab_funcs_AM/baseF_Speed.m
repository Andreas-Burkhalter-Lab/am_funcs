function [F08,F0]=baseF_Speed(Fraw,F,rotation,forward,varargin)
if nargin>4
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>5
        env=varargin{3};
    else env='';
    end
else
    saveDir='F:\MA Data\Interneurons\';
    MouseID='This Mouse';
    
end

F0=zeros(size(Fraw,2),1);
F08=zeros(size(Fraw,2),1);
corrs=F0;
numframes=length(Fraw);
numwindow=min(500,numframes/50);

for j=1:size(Fraw,2)
    junk=Fraw(:,j);
    
    window=round(numframes/numwindow);
    junk2=zeros(size(junk));
    for k=1:length(junk)
        cut=junk(max(1,k-window):min(numframes,k+window));
        cutsort=sort(cut);
        a=round(length(cut)*.08);
        junk2(k)=cutsort(a);
    end
    F08(j)=mean(junk2);
    if ~isnan(nanmean(Fraw(:,j)))
    sm=smooth(Fraw(:,j),10);
    else 
        sm=nan;
    end
    F0(j)=min(sm);
    ctemp=corrcoef(F(:,j),sqrt(forward(:,j).^2+rotation(:,j).^2));
    corrs(j)=ctemp(2,1);
end

figure
subplot(1,2,1)
scatter(F0,corrs)
title('F0 vs correlation with speed')
xlabel('F0')
subplot(1,2,2)
scatter(F08,corrs)
title('F08 vs correlation with speed')
xlabel('F08')



saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'F0 vs correlation with speed', '.jpg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'F0 vs correlation with speed', '.fig']);
