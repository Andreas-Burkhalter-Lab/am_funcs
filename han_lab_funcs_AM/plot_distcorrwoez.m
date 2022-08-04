function plot_distcorr(mouse,name,varargin)
if nargin>2
    mouseNums=varargin{1};
    expnum=varargin{2};
%     Fs=varargin{3};
    if nargin>4
        saveDir=varargin{3};
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
distsall=[];
corrsall=[];
figure
hold on
leg{1}=[];
ii=1;
for m=mouseNums
    for e=expnum{m}
        if e<=length(mouse(m).Dists)
        d=mouse(m).Distswoez{e};
        c=mouse(m).Corrswoez{e};
        distsall=[d,distsall];
        corrsall=[c,corrsall];
        if ~isempty(d)
        scatter(d,c);

        leg{ii}=['Mouse ',num2str(m),' exp ',num2str(e)];
        ii=ii+1;
        end
        end
    end
    
end
legend(leg)
% scatter((distsall),corrsall)
title('Dist vs correlation wo EZ')
xlabel('Distance')
ylabel('Correlation')

[saveDir,'\', 'Dist vs correlation woEZ',' All Cells',name,'.jpg']
        saveas(gca,[saveDir,'\', 'Dist vs correlation woEZ',' All Cells',name,'.jpg']);
        savefig([saveDir,'\', 'Dist vs correlation woEZ', ' All Cells ',name,'.fig']);
        