function [noise,IN,PN]=split_cells_CPP(Fc,varargin)
if nargin>1
    saveDir=varargin{1};
    MouseID=varargin{2};
else
    saveDir='H:\MA Data\Mixed\PubFigs';
    MouseID='This guy';
end
noise_thresh=.5;
IN_thresh=2;


if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
[noise,IN,PN]=plot_check(Fc,[],['Fc File',' Split Cells'],saveDir,MouseID);
% Fc(:,noise)=nan(length(Fc),sum(noise));

% close all
    function [noise,IN,PN]=plot_check(F,~,name,saveDir,MouseID)
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
%         ts=tight_subplot(1,2,[.01 .01],[.01 .01],[.01,.01]);
%         s1=subplot(size(F,2),2,1:2:(2*size(F,2)));
%         axes(ts(1));
        plot(bsxfun(@plus,F,1:size(F,2)))
        set(gca,'Ytick',1:size(F,2),'YTickLabel',strsplit(num2str((skewness(F))),' '))
%         for cell=1:size(F,2)
%            ss(cell)=subplot(size(F,2),2,(2*size(F,2)+2)-2*cell);
%            histogram(F(:,cell))
%            set(gca,'view',[90,-90])
%            xlabel(num2str(skewness(F(:,cell))),'Rotation',-90)
% %            xlabel('test')
%            set(get(gca,'XLabel'),'Rotation',0)
%         end
%         linkaxes([ts(1),ts(2)],'y')
%         legend('p2r','p2p','std','var','rms','rssq')
        suptitle(name)
        if ~exist([saveDir,'\','Cell Tests' ,'\'],'dir')
            mkdir([saveDir,'\','Cell Tests' ,'\']);
        end
        saveas(gca,[saveDir,'\','Cell Tests' ,'\', MouseID, 'Cell Test', name , '.jpg']);
        savefig([saveDir,'\','Cell Tests','\', MouseID, 'Cell Test', name ,'.fig']);
        PN=find(skewness(F)>IN_thresh);
        IN=find(skewness(F)<.6 & skewness(F)>noise_thresh);
        noise=find(skewness(F)<noise_thresh);
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        subplot(1,3,1)
        plot(bsxfun(@plus,F(:,PN),1:size(F(:,PN),2)))
        set(gca,'Ytick',1:size(F(:,PN),2),'YTickLabel',strsplit(num2str(PN)))
        title('Pyramidal')
        subplot(1,3,2)
        plot(bsxfun(@plus,F(:,IN),1:size(F(:,IN),2)))
        set(gca,'Ytick',1:size(F(:,IN),2),'YTickLabel',strsplit(num2str(IN)))
        title('Interneuron')
        subplot(1,3,3)
        plot(bsxfun(@plus,F(:,noise),1:size(F(:,noise),2)))
        set(gca,'Ytick',1:size(F(:,noise),2),'YTickLabel',strsplit(num2str(noise)))
        title('Noise')
        saveas(gca,[saveDir,'\','Cell Tests' ,'\', MouseID, 'Split', name , '.jpg']);
        savefig([saveDir,'\','Cell Tests','\', MouseID, 'Split', name ,'.fig']);

        
    end
end