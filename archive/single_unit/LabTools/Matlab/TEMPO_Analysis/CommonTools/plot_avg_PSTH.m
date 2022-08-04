%This function loads the file 'filename' which is a cumulative record of
%spike times.  The function returns a plot of the cumulative PSTH.  CMA
%06/28/06.
%function plot_avg_PSTH(filename)
spikes = ALSpikes_OptimizedSF;
[a,b] = size(spikes);
spikes = reshape(spikes, 1, (a*b));
figure
hist(spikes,60)
%histogram has 50 ms bins
title('Optimized SF PSTH - AL')
xlabel('Time (s)'), ylabel('Spikes')