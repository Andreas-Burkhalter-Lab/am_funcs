%--------------------------------------------------------------------------------------------------------
% SpikeBinner.m: Tool for binning spike data into larger bins.  This function assumes that spike bins have
%   values > 0 if at least one spike occurred, and is intended mainly for use in binning spikes into
%   histograms.  The bin widths are
%   assumed to be in milliseconds.  Offset_bin can be used to have the output 'bins' offset by a certain
%   amount of time. 
%   Typical usage: [bins, counts] = SpikeBinner(spike_data, 1, 20, 200)
%       This will take spike data sampled at 1ms resolution and bin the spikes into 20ms bins,
%       with the output bin times starting at -200ms.  The binned histogram can easily be plotted using:
%       bar(bins, counts).  NOTE: raw_event_stream should be a COLUMN VECTOR, counts is a COLUMN VECTOR
%   GCD, starting 6/20/01

function [bins, counts] = SpikeBinner(raw_event_stream, raw_bin_width, resampled_bin_width, offset_bin)

    %convert to column vector if not one already
    if ( (size(raw_event_stream,2)>1) & (size(raw_event_stream,1)==1) )	% a row vector, so transpose
        raw_event_stream = raw_event_stream';
    end
    %first, find the maximum number of spikes in any original bin
    %the above code can be shortened... - BJP
    max_spikes_per_bin = max(raw_event_stream);
    
    %now, accumulate the event times for all spikes, counting times multiply if the original bin
    %value exceeds 1
    spike_times = [];
    for j=1:max_spikes_per_bin
        st = find(raw_event_stream == j);
        temp = 1;
        while (temp <=j)
            spike_times = [spike_times; st];
            temp = temp + 1;
        end        
    end

    %now, add in the offset and multiply by the raw bin width (usually 1)
    spike_times = (spike_times - offset_bin).*raw_bin_width;
    
    num_raw_bins = length(raw_event_stream);
    num_resampled_bins = floor(num_raw_bins*raw_bin_width/resampled_bin_width);
    
    begin_time = (0 - offset_bin)*raw_bin_width;
    end_time = (length(raw_event_stream) - offset_bin)*raw_bin_width;
    edges = begin_time:resampled_bin_width:end_time;  
    bins = edges  +  resampled_bin_width/2;
    
    counts = [];
    if ~(isempty(spike_times) )
        [counts] = histc(spike_times, edges);
    end
    if (isempty(counts) )
        counts = 0 * bins';   % make sure vector is full of zeros and not empty for output
    end   
    
    
    %convert to column vector if not one already
    if ( (size(counts,2)>1) & (size(counts,1)==1) )	% a row vector, so transpose
        counts = counts';
    end
  
    % bins and counts appear to have one more index than needed
    %also, keep vector orientation for output same as input
    % histc returns a vector whose last index contains the number of counts matching the last edge.
    %   we deleted that
    counts = counts(1:end-1);
    % bins contains indices for time center of each bin
    bins = bins(1:end-1)';
    
return;