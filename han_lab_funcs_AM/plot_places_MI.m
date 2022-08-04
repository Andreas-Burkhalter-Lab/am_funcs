%Plots feature_count place fields with MI in x
function plot_places_MI(feature_count, MI,saveDir,MouseID,label,varargin)
if nargin>5
    placeRow=varargin{1};
end
if ~exist([saveDir,'\',MouseID,'\'],'dir')
    mkdir([saveDir,'\',MouseID,'\'])
end
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
    MIinds=strsplit(num2str(MI));

    set(gca, 'XTick', 1:size(feature_count,2), 'XTickLabel', MIinds, 'XTickLabelRotation', 90);
    subplot(2,1,2)
    feature_sorted=feature_count(:,inds);
    imagesc(feature_sorted);
    set(gca,'Ydir','normal')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID,' Places ',label,' .jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Places ',label,' .fig']);
    
    sinds=strsplit(num2str(inds));
    set(gca, 'XTick', 1:size(feature_count,2), 'XTickLabel', sinds, 'XTickLabelRotation', 90);
    if exist('placeRow','var')
        import mlreportgen.dom.*;
        append(placeRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID,' Places ',label,' .jpg'])))
    end
        
end