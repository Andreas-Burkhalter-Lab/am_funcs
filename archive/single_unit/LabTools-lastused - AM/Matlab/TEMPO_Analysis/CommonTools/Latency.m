% function [lt] = Latency(spikes,stim_time)
% 
% Latency.m:  Computes response latency from previously-aligned spike trains using
% three methods (described below) and returns them in ms as row vector 'lt'.
% Raw data 'spikes' matrix should be 3 dimensional: SpikeChan x Time x
% Trial, just like in the data.spike_data structure.  You can send the
% entire data.spike_data structure, but you might want to pass only those
% conditions in which you expect a strong response.
% The 'stim_time' is the stimulus onset time in ms.  You must identify this time and align
% your spike trains to it, but be sure to leave 100-250ms before stimulus
% onset to use as a baseline.  (If you have no variable delays, your
% data.spike_data may already be aligned.)
%
% Method One finds the first bin to hit 50% of the peak response in a mass histogram of responses.
% (Smith 2005 in MT)
%
% Method Two fits all the time preceding stimulus onset with a Poisson, and then finds the first of 
% three consecutive bins where the count has a p<0.01 chance of coming from that distribution.
% (Maunsell 1992 in V1 and Bisley 2004 in PPC)
%
% Method Three calculates a mean and stddev firing rate from the time preceding stimulus
% onset and then finds the first bin in the mass histogram that exceeds 2*sd.
%
% JWN 06/09/05
% JWN 06/13/05

function [lt] = Latency(spikes,stim_time)

% Initialize with zeros
lt = [0 0 0];

% Both methods require PSTHs
% Using 4ms bins gives you many noisy bins vs. fewer bins with better STN
bin_width = 4;
[bins, counts] = SpikeBinner(sum(spikes,3),1,bin_width,stim_time);  % Bins and shifts data so stim_time is 0
            
% Method One
peak = max(counts);  % Find the peak response
tmp = find(counts>max(counts)/2 & bins>0);  % Find all bins with >50% peak after stim_time
lt(1) = bins(tmp(1));  % Find the time of the first >50% bin after stim_time

% Method Two
lambda = poissfit(counts(find(bins<0)));  % Fit preceding time baseline to a Poisson
% Compute the Poisson cumulative distribution function with parameter
% lambda at the values in the bins following stim_time
tmp = poisscdf(counts,lambda)>=.99;  % Thresholded at p<0.01
tmp2 = find(conv(tmp,[1 1 1])==3)-2;  % Find sequences of three consecutive qualifying bins
lt(2) = bins(tmp2(min(find(bins(tmp2)>0))));  % Find the time of the first bin from the above set

% Method Three
baseline_mean = mean(counts(find(bins<0)));  
baseline_stddev = std(counts(find(bins<0)));
tmp = find(counts>(baseline_mean+3*baseline_stddev) & bins>0);  % Find the bins that exceeds mean+2*sd after stim_time
lt(3) = bins(tmp(1));
return;