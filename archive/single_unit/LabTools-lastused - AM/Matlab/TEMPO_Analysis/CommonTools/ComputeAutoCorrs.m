%----------------------------------------------------------------------------------------------------------
%-- ComputeAutoCorrs.m: This function computes the autocorrelation for each spike train and each trial. The
%-- range of lags is determined by -MaxLag->MaxLag. The autocorrelation data are returned for each trial and 
%-- each spike channel in autocorr_matrix(trial#,lag,chan#).  The lag values are also returned in lags.
%-- GCD, 6/5/01
%----------------------------------------------------------------------------------------------------------
function [autocorr_matrix, lags] = ComputeAutoCorrs(data, n_trials, start_code, stop_code, StartOffset, StopOffset, MaxLag);

TEMPO_Defs;	%some defines that we'll need

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{SPIKE_DB};	%for convenience
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

%Here, I extract spikes over the period from start_index to stop_index for every trial
%In many cases, start_index and stop_index would be the same for each trial, and the loops would not be needed
%But, to be general, I am doing it this way so that each trial could have different start- and stop_ indices

%instead of using the number of channels in the header, switch to calculating the number of channels from the size of the data array
%-JDN 11/29/00
size_data = size(data.spike_data);
num_chan = size_data(1);

%preallocate the autocorr_matrix for speed
autocorr_matrix = ones(n_trials, 2*MaxLag+1, num_chan);

for j = 1:num_chan		%for each spike channel
    for i = 1:n_trials		%for each trial        
        start_eventbin = find(data.event_data(1,:,i) == start_code);
        stop_eventbin = find(data.event_data(1,:,i) == stop_code);
        
        % convert to spike bins, add offsets (already in spike bin format)
        start_spikebin = floor (start_eventbin*(event_bin_width/spike_bin_width)) + StartOffset;
        stop_spikebin = floor (stop_eventbin*(event_bin_width/spike_bin_width)) + StopOffset;
        
        spikes = data.spike_data(j,start_spikebin:stop_spikebin,i);
        [counts, lags] = xcorr(spikes, spikes, MaxLag);
        
        autocorr_matrix(i,:,j) = counts;    
    end
end

return;
