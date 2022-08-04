function check_cells(Fc,varargin)
if nargin>1
    saveDir=varargin{1};
    MouseID=varargin{2};
else
    saveDir=pwd;
    MouseID='This guy';
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
demF=bsxfun(@minus,Fc,mean(Fc));
normF=demF/(max(demF(:)));
% 
% Fcell=mat2cell(F,size(F,1),ones(1,size(F,2)));
% % dFcell=mat2cell(demF,size(F,1),ones(1,size(F,2)));
% nFcell=mat2cell(normF,size(F,1),ones(1,size(F,2)));
% Fccell=mat2cell(Fc,size(F,1),ones(1,size(F,2)));
% 
% plot_check(normF,Fcell,'F File',saveDir,MouseID)
% plot_check(normF,dFcell,'Demeaned F File',saveDir,MouseID)
plot_check(normF,[],'Normalized F File',saveDir,MouseID)
plot_check(Fc,[],'Fc File',saveDir,MouseID)

    function plot_check(F,~,name,saveDir,MouseID)
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         ts=tight_subplot(1,2,[.01 .01],[.01 .01],[.01,.01]);
        s1=subplot(size(F,2),2,1:2:(2*size(F,2)));
%         axes(ts(1));
        plot(bsxfun(@plus,F,1:size(F,2)))
        for cell=1:size(F,2)
           ss(cell)=subplot(size(F,2),2,(2*size(F,2)+2)-2*cell);
           histogram(F(:,cell))
           set(gca,'view',[90,-90])
           xlabel(num2str(skewness(F(:,cell))))
        end
%         linkaxes([ts(1),ts(2)],'y')
%         legend('p2r','p2p','std','var','rms','rssq')
        suptitle(name)
        if ~exist([saveDir,'\','Cell Tests' ,'\'],'dir')
            mkdir([saveDir,'\','Cell Tests' ,'\']);
        end
        saveas(gca,[saveDir,'\','Cell Tests' ,'\', MouseID, 'Cell Test', name , '.jpg']);
        savefig([saveDir,'\','Cell Tests','\', MouseID, 'Cell Test', name ,'.fig']);
    end
end