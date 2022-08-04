function [out]=mean_and_sem(data)

m=nanmean(data);
sem=nanstd(data)/sqrt(length(data));
out=[m sem];