% ytrials: param values correspond to rows, repetitions correspond to
% columns

% set breakpoint in [Parameter]TuningCurve.m after y_interp is calculated

clear ytrials match

% turn on if parameter is direction or orientation tuning
dir_orient_tuning = 1;

switch dir_orient_tuning
    case 0
        ytrials = zeros(length(px),round(length(plot_y) / [length(px)-1]));
        for i = 1:length(px)-1
            match = plot_x==px(i);
            ytrials(i,:) = plot_y(match);
        end
    case 1
        ytrials = zeros(length(px),round(length(plot_y) / [length(px)]));
        for i = 1:length(px)
            match = plot_x==px(i);
            ytrials(i,:) = plot_y(match);
        end
end

save('16050 c3 r13 of coherence tuning','py','px','plot_x','plot_y','perr','ytrials','x_interp','y_interp')

