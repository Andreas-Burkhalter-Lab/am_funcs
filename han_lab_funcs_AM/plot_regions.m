%Plots feature_count place fields with MI in x
function plot_regions(feature_count,saveDir,MouseID,label,varargin)
if nargin>4
    placeRow=varargin{1};
end
if ~exist([saveDir,'\',MouseID,'\'],'dir')
    mkdir([saveDir,'\',MouseID,'\'])
end
    feature_count=bsxfun(@rdivide,feature_count,max(feature_count,[],1));
    [~,maxind] = max(feature_count);
    [~, inds] = sort(maxind);
    %Plot location of max response for each cell
    figure();
    stem(maxind);
    %Plot tuning curves and sorted tuning curve
    f1=figure();
    set(f1,'Position', [5 15 1900 980]);
    subplot(2,1,1)
    imagesc(feature_count)
    set(gca,'Ydir','normal')
%     MIinds=strsplit(num2str(MI));

    set(gca, 'XTick', 1:size(feature_count,2));
    subplot(2,1,2)
    feature_sorted=feature_count(:,inds);
    imagesc(feature_sorted);
    set(gca,'Ydir','normal')
    jpgName=[saveDir,'\',MouseID,'\', MouseID,' Regions ',label,' .jpg'];
    saveas(gca,jpgName);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Regions ',label,' .fig']);
    
    set(gca, 'XTick', 1:size(feature_count,2));
    if exist('placeRow','var')
        import mlreportgen.dom.*;
        append(placeRow,TableEntry(Image(jpgName)))
    end
        
end