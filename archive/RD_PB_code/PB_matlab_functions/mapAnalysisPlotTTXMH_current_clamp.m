function mapAnalysisPlotTTXMH_current_clamp(map)

%%% AM edited to use with current clamp

% makes a figure with plots with the analyzed parameters of a map.
% LTP 2006
[r c]=size(map.pattern{1});
dx=map.xSpacing;
dy=map.ySpacing;
xdata = [0:dx:(c-1)*dx] - (c-1)*dx/2;
ydata = [0:dy:(r-1)*dy] - (r-1)*dy/2;; 
f=figure('Position',[ 884    13   710   698]);

X=map.soma1Coordinates(1); %-map.xPatternOffset; MH attempt to correct soma plotting
Y=map.soma1Coordinates(2); %-map.yPatternOffset; MH attempt to correct soma plotting


subplot(3,3,1)
title('mean');
hax=gca;
h=imagesc(map.analysis.mean);

colorbar
set(h, 'XData', xdata, 'YData', ydata);
set(hax, 'XLim', [min(xdata)-dx/2 max(xdata)+dx/2], 'YLim', [min(ydata)-dy/2 max(ydata)+dy/2])
hold on
hplt = plot(X, -Y, 'k^');
set(hplt,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
title('mean');
hold off

subplot(3,3,2)
title('maximum');       %% AM 4/7/15 changed minimum to maximum
hax=gca;
h=imagesc(map.analysis.maximum);       %% AM 4/7/15 changed map.analysis.minimum to map.analysis.maximum

colorbar
set(h, 'XData', xdata, 'YData', ydata);
set(hax, 'XLim', [min(xdata)-dx/2 max(xdata)+dx/2], 'YLim', [min(ydata)-dy/2 max(ydata)+dy/2])
title('maximum');       %% AM 4/7/15 changed minimum to maximum
hold on
hplt = plot(X, -Y, 'k^');
set(hplt,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
hold off

subplot(3,3,3)

hax=gca;
h=imagesc(map.analysis.onset);

colorbar
set(h, 'XData', xdata, 'YData', ydata);
set(hax, 'XLim', [min(xdata)-dx/2 max(xdata)+dx/2], 'YLim', [min(ydata)-dy/2 max(ydata)+dy/2])
title('onset');
hold on
hplt = plot(X, -Y, 'k^');
set(hplt,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
hold off

subplot(3,3,4)

hax=gca;
h=imagesc(map.analysis.riseTime1090);

colorbar
set(h, 'XData', xdata, 'YData', ydata);
set(hax, 'XLim', [min(xdata)-dx/2 max(xdata)+dx/2], 'YLim', [min(ydata)-dy/2 max(ydata)+dy/2])
title('rise time');
hold on
hplt = plot(X, -Y, 'k^');
set(hplt,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
hold off


%%% AM commented out 4/7/15
% % % % % % % % % % subplot(3,3,5)
% % % % % % % % % % 
% % % % % % % % % % hax=gca;            
% % % % % % % % % % h=imagesc(map.analysis.integral);
% % % % % % % % % % %colormap(flipud(jet2))
% % % % % % % % % % colorbar
% % % % % % % % % % set(h, 'XData', xdata, 'YData', ydata);
% % % % % % % % % % set(hax, 'XLim', [min(xdata)-dx/2 max(xdata)+dx/2], 'YLim', [min(ydata)-dy/2 max(ydata)+dy/2])
% % % % % % % % % % title('charge');
% % % % % % % % % % hold on
% % % % % % % % % % hplt = plot(X, -Y, 'k^');
% % % % % % % % % % set(hplt,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
% % % % % % % % % % hold off


subplot(3,3,7)

hist(map.analysis.onset(:),30);
title('onset');

subplot(3,3,8)

hist(map.analysis.riseTime1090(:),30);
title('rise time');

subplot(3,3,9)
ax=gca;
set(ax,'Visible','Off');

titleStr = [map.experimentNumber];
text('String', titleStr, 'Units', 'Normalized', 'Position', [0.33 0.17], ...
    'FontSize', 12, 'FontWeight', 'Bold');

plotMappingArrayOverImage(map)

figure()
subplot(1,2,1)
title('mean');
hax=gca;
h=imagesc(map.analysis.mean);
colorbar
set(h, 'XData', xdata, 'YData', ydata);
set(hax, 'XLim', [min(xdata)-dx/2 max(xdata)+dx/2], 'YLim', [min(ydata)-dy/2 max(ydata)+dy/2])
hold on
hplt = plot(X, -Y, 'k^');
set(hplt,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
title('mean');
hold off
xlabel('Distance (\mum)')
ylabel('Distance (\mum)')

%%% AM commented out 4/7/15
% % % % % subplot(1,2,2)
% % % % % hax=gca;
% % % % % h=imagesc(map.analysis.integral);
% % % % % %colormap(flipud(jet2))
% % % % % colorbar
% % % % % set(h, 'XData', xdata, 'YData', ydata);
% % % % % set(hax, 'XLim', [min(xdata)-dx/2 max(xdata)+dx/2], 'YLim', [min(ydata)-dy/2 max(ydata)+dy/2])
% % % % % title('charge');
% % % % % hold on
% % % % % hplt = plot(X, -Y, 'k^');
% % % % % set(hplt,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
% % % % % hold off
% % % % % xlabel('Distance (\mum)')
% % % % % ylabel('Distance (\mum)')
end