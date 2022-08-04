%%% test bandwidth values for violin plot
% updated 2020/9/21 on thermaltake



bws = logspace(0.0005,0.06,20); % values for pakan loc index

ops.newfig = 0;
ops.close_all = 0;
for i = 1:20
    subplot(4,5,i)
    bw = log10(bws(i));
    ops.violin_bandwidth = bw;
    tuning_stats
    title(num2str(bw))
end