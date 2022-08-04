function plot_places(feature_count, varargin)
%sort by point of max response
%     feature_count(:,[7 53 105])=[];
if nargin>2
    prior=varargin{2};
else prior=[];
end
if nargin>3
    pname=varargin{3};
    name=varargin{4};
    type=varargin{5};
end
if nargin>6
    part=varargin{6};
    import mlreportgen.dom.*;
end
[~,maxind] = max(feature_count);
[maxes, inds] = sort(maxind);

%Plot tuning curves and sorted tuning curve
f1=figure();
set(f1,'Position', [5 15 1900 980]);
subplot(2,1,1)
imagesc(feature_count)
%     colormap('hot')
title('Cell Location Preference (Mean normalized dF at each loc)')
colorbar()
ylabel('Position')
set(gca,'YTick',[12 24],'YTickLabel',{'A','B'},'ydir','normal')
subplot(2,1,2)
feature_sorted=feature_count(:,inds);
imagesc(feature_sorted);
title('Sorted Location Preference')
ylabel('Position')
xlabel('Cells')
set(gca,'YTick',[12 24],'YTickLabel',{'A','B'},'ydir','normal')

colormap('jet')
colorbar()

%     colormap('hot')
sinds=strsplit(num2str(inds));
set(gca, 'XTick', 1:size(feature_count,2), 'XTickLabel', sinds, 'XTickLabelRotation', 90);
if nargin>2
    saveas(gca,[pname,name,'\',name,' PlaceSpecificity', type,'.jpg']);
        saveas(gca,[pname,name,'\',name,' PlaceSpecificity', type,'.svg']);

end

if nargin>4
    I=Image([pname,name,'\',name,' PlaceSpecificity', type,'.jpg']);
    append(part,TableEntry(I));
end
%Plot location of max response for each cell
figure('units','normalized', 'Position', [.01 .05 .95 .86]);
subplot(1,3,1)
bar(prior)
title('Mouse Location Residency Time')
set(gca,'XTick',[12 24],'XTickLabel',{'A','B'})
ylabel('Proportion Time')
xlabel('Region')
subplot(1,3,2)
histogram(maxind)
title('Cell Location Preference')
xlabel('Binned position (5cm)')
ylabel('Number of Cells Preffering')
set(gca,'XTick',[12 24],'XTickLabel',{'A','B'})
subplot(1,3,3)
bar(([sum(maxind>1 & maxind<12) sum(maxind>12 & maxind<24) sum(maxind>24 & maxind<37)]/length(maxind))')
title('Cell Location Preference Regions')
ylabel('Proportion Time')
xlabel('Region')
set(gca,'XTick',1:3,'XTickLabel',{'A','Mid','B'})
if nargin>2
    saveas(gca,[pname,name,'\',name,' HistogramPlaceSpecificity', type,'.jpg']);
        saveas(gca,[pname,name,'\',name,' HistogramPlaceSpecificity', type,'.svg']);

end
if nargin>4
    I=Image([pname,name,'\',name,' HistogramPlaceSpecificity', type,'.jpg']);
    append(part,TableEntry(I));
end
if nargin>1 && ~isempty(varargin{1})
    feature_count_2=varargin{1};
    %sort by point of max response
    [~,maxind_2] = max(feature_count_2);
    [maxes_2, inds_2] = sort(maxind_2);
    %Plot location of max response for each cell
    figure();
    subplot(1,2,1)
    stem(maxind_2);
    subplot(1,2,2)
    histogram(maxind_2)
    
    %Plot tuning curves and sorted tuning curve
    f1=figure();
    set(f1,'Position', [5 15 1900 980]);
    subplot(2,1,1)
    imagesc(feature_count_2)
    colormap('jet')
    subplot(2,1,2)
    feature_sorted_2=feature_count_2(:,inds);
    imagesc(feature_sorted_2);
    sinds_2=strsplit(num2str(inds));
    set(gca, 'XTick', 1:size(feature_count_2,2), 'XTickLabel', sinds_2, 'XTickLabelRotation', 90);
end
end



%
%
% plot_places(feature_count2);
% h=gca;
% saveas(h,[saveDir,'\',MouseID,'\', MouseID,' Test Normalized by Prior and Gaussian Activity Map.jpg']);
%
%
% plot_places(feature_count1);
% h=gca;
% saveas(h,[saveDir,'\',MouseID,'\', MouseID,' Training Normalized by Prior and Gaussian Activity Map.jpg']);
%
%
% plot_places(feature_countb1);
% h=gca;
% saveas(h,[saveDir,'\',MouseID,'\', MouseID,' Training Binary Not Normalized Activity Map.jpg']);
%
%
% plot_places(feature_countb2);
% h=gca;
% saveas(h,[saveDir,'\',MouseID,'\', MouseID,' Test Binary Not Normalized Activity Map.jpg']);
%
% close all;