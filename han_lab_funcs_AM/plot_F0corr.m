function plot_F0corr(mouse,varargin)
if nargin>1
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>4
        saveDir=varargin{4};
    else     saveDir='F:\MA Data\Interneurons\PubFigs';
    end
    
else
    mouseNums=1:length(mouse);
    for m=mouseNums
        expnum{m}=1:length(mouse(m).Falls);
    end
    Fs=4*ones(size(mouseNums));
    saveDir='F:\MA Data\Interneurons\PubFigs';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
F0all=[];
F08all=[];
corrsall=[];
for m=mouseNums
    for e=expnum{m}
        
        F1=mouse(m).Falls{e};
        f1=mouse(m).Forwards{e};
        r1=mouse(m).Rotations{e};
        F0=mouse(m).F0{e};
        F08=mouse(m).F08{e};
        corrs=zeros(size(F0,1),1);
        for j=1:size(F0)
            ctemp=corrcoef(F1(:,j),sqrt(f1(:,j).^2+r1(:,j).^2));
            corrs(j)=ctemp(2,1);
        end
        F0all=[F0;F0all];
        F08all=[F08;F08all];
        corrsall=[corrs;corrsall];
        m
        e
        scatter(F0,corrs)
    end
    
end
figure
subplot(1,2,1)
scatter(F0all,corrsall)
title('F0 vs correlation with speed')
xlabel('F0')
subplot(1,2,2)
scatter(F08all,corrsall)
title('F08 vs correlation with speed')
xlabel('F08')

        saveas(gca,[saveDir,'\', 'F0 speed corr',num2str(m),',',num2str(e),' All Cells','.jpg']);
        savefig([saveDir,'\', 'F0 speed corr',num2str(m),',',num2str(e), ' All Cells ','.fig']);
        