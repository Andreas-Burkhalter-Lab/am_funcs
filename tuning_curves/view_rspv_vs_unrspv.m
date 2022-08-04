%%% view responsive vs unresponsive cells' dff
%%% first load dff_ file and corresponding tuning struct output from get_tuning_curves.m
roirangeind = 1;
rows = 7;
cols = 4;

tuning = p.tuning; % select tuning table
scope_events = res.scope_events;

% find responsive vs unresponsive rois
is_rspv = tuning.sf_sgnf | tuning.tf_sgnf | tuning.orient_sgnf;
rspv = find(is_rspv);
unrspv = find(~is_rspv);

nplots = rows*cols;
roirange = [1:nplots] + nplots*[roirangeind-1];
close all


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:nplots
    innerind = roirange(i);
    if i>length(rspv)
        break
    end
    tuningrowind = rspv(innerind);
    subplot(rows,cols,i)
    f_to_plot = res.scope_events.dff(:,tuningrowind);
    hold on
    plot(res.scope_events.stim_present*0.5*max(f_to_plot)) % plot stim timing
    plot(f_to_plot)
    hold off
    title(['tuning row ' num2str(tuningrowind)])
end
suptitle('responsive')

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:nplots
    innerind = roirange(i);
    if innerind>length(unrspv)
        break
    end
    tuningrowind = unrspv(innerind);
    subplot(rows,cols,i)
    f_to_plot = res.scope_events.dff(:,tuningrowind);
    hold on
    plot(res.scope_events.stim_present*0.5*max(f_to_plot)) % plot stim timing
    plot(f_to_plot)
    hold off
    title(['tuning row ' num2str(tuningrowind)])
end
suptitle('unresponsive')