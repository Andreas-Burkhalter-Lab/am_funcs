function [out]=mean_and_std(data)

m=mean(data);
sem=std(data);
out=[m, sem];