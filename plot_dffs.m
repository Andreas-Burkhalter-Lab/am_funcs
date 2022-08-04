roirangeind = 2;
rows = 7;
cols = 4;


plot_dff_precut = 0; 
plot_dff_postcut = 1;



nplots = rows*cols;
roirange = [1:nplots] + nplots*[roirangeind-1];
close all

if plot_dff_precut
    figure('units','normalized','outerposition',[0 0 1 1])
    for roirange(i) = 1:nplots
        if i>size(dff,2)
            break
        end
        subplot(rows,cols,i)
        plot(dff(:,roirange(i)))
        title(num2str(roirange(i)))
    end
    suptitle('dff')
end


%%% 
if plot_dff_postcut
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:nplots
        if roirange(i)>size(res.scope_events.dff,2)
            break
        end
        subplot(rows,cols,i)
        plot(res.scope_events.dff(:,roirange(i)))
        title(num2str(roirange(i)))
    end
    suptitle('dff with bad frames cut')
end