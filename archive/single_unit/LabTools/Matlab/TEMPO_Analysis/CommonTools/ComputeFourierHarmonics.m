%----------------------------------------------------------------------------------------------------------
%-- ComputeFourierHarmonics.m: This function computes the Fourier responsee harmonic for each spike channel and each
%--	trial, over a period of time specified by start_code and stop_code, which may occur at different
%--	times on different trials.  Code returns the amplitude of the first 10 harmonics. GCD, 7/29/02
%----------------------------------------------------------------------------------------------------------
function harmonics = ComputeFourierHarmonics(all_data, n_trials, start_code, stop_code, StartOffset, StopOffset);

TEMPO_Defs;	%some defines that we'll need

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = all_data.htb_header{SPIKE_DB};	%for convenience
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed)

h = all_data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

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
        
        tt = ((start_spikebin:stop_spikebin) - start_spikebin).*spike_bin_width;
        raster = all_data.spike_data(j,start_spikebin:stop_spikebin,i);
        DC_subtract = 1;
        [freq, ampl] = FourierTransform_1D(tt, raster, length(raster), DC_subtract, 0);
        for k = 1:10
            harmonics(j,i,k) = ampl(k+1);
        end
        %harmonics has dimensions (channel#, trial#, harmonic#)
    end
end

return;
