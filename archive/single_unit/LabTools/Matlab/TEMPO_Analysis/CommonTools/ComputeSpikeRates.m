%----------------------------------------------------------------------------------------------------------
%-- ComputeSpikeRates.m: This function computes the average firing rate for each spike channel and each
%--	trial, over a period of time specified by start_code and stop_code, which may occur at different
%--	times on different trials.  GCD, 12/27/99
%----------------------------------------------------------------------------------------------------------
function spike_rates = ComputeSpikeRates(all_data, n_trials, start_code, stop_code, StartOffset, StopOffset);

TEMPO_Defs;	%some defines that we'll need


%first, we'll need the bin_width (in sec) of our spike raster + events log
h = all_data.htb_header{SPIKE_DB};	%for convenience
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

h = all_data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

%Here, I count spikes over the period from start_index to stop_index for every trial
%In many cases, start_index and stop_index would be the same for each trial, and the loops would not be needed
%But, to be general, I am doing it this way so that each trial could have different start- and stop_ indices

%instead of using the number of channels in the header, switch to calculating the number of channels from the size of the data array
%-JDN 11/29/00
size_data = size(all_data.spike_data);
num_chan = size_data(1);
for j = 1:num_chan		%for each spike channel
    for i = 1:n_trials		%for each trial        
        start_eventbin = find(all_data.event_data(1,:,i) == start_code);
        stop_eventbin = find(all_data.event_data(1,:,i) == stop_code);
        
        % convert to spike bins, add offsets (already in spike bin format)
        start_spikebin = floor (start_eventbin*(event_bin_width/spike_bin_width)) + StartOffset;
        stop_spikebin = floor (stop_eventbin*(event_bin_width/spike_bin_width)) + StopOffset;
        
        count = sum(all_data.spike_data(j,start_spikebin:stop_spikebin,i), 2);
        spike_rates(j,i) = count / ((stop_spikebin - start_spikebin)*spike_bin_width);
        %spike_counts has dimensions (channel#, trial#), and each value is a spike count over the range of sample indices given by start_index, stop_index
    end
end

return;
