%----------------------------------------------------------------------------------------------------------
%-- CheckTimeOffset.m: This function checks the offset times entered into the Tempo GUI to ensure that they
%--	are within the time ranges of the database and outputs (in spike bins, offsets good for
%-- 	all trials.  Output the bin number corresponding to the specified event and time offset for each trial
%--   I also have it output the actual bins (without offset) for each trial - BJP, 5/1/01
%----------------------------------------------------------------------------------------------------------
function [StartOffsetIndex, StopOffsetIndex, start_eventbin, stop_eventbin] = CheckTimeOffset(all_data, num_trials, start_code, stop_code, StartOffsetTime, StopOffsetTime, UseSyncPulses);

%index for last original spike channel - usually sync pulses
global num_recorded_spike_channels;

TEMPO_Defs;	%some defines that we'll need

h = all_data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

StartOffsetIndex = (StartOffsetTime/1000)/event_bin_width; % convert times to ms then get indices
StopOffsetIndex = (StopOffsetTime/1000) / event_bin_width;

%Checks start and stop offsets againsts restraints of each trial in even stream

listtext = [];
for i = 1:num_trials		%for each trial
    start_eventbin(i) = find(all_data.event_data(1,:,i) == start_code);   
    stop_eventbin(i) = find(all_data.event_data(1,:,i) == stop_code);    
    
%    if (~isempty(all_data.spike_data) ) & (sum(all_data.spike_data(size(all_data.spike_data, 1),:,i) > 0) & UseSyncPulses)
    if (~isempty(all_data.spike_data) ) & (sum(all_data.spike_data(num_recorded_spike_channels ,:,i) > 0) & UseSyncPulses)
        %make sure sync pulses are not empty.
        if (start_code == VSTIM_ON_CD & size(all_data.spike_data, 1) > 1)
            %if startcode is stim start and have sync pulses, overwrite tempo drop code times
            %check for actual starting and ending frames of visual stimulus
            %sync pulses are stored in last channel of spike_data
    	    start_eventbin(i) = min( find(all_data.spike_data(num_recorded_spike_channels,:,i) == 1) );
            %% next lines check for spurious ones in the sync channel,
            %% present in m5c118r5 data set - BJP 2/18/03
            if start_eventbin(i) < 500
                
               sync_times = find(all_data.spike_data(size(all_data.spike_data, 1),:,i) == 1);
               start_eventbin(i) = min(sync_times(sync_times > 500) );
            end       
        end
    
        if (start_code == VSTIM_OFF_CD & size(all_data.spike_data, 1) > 1)
            %if startcode is stim start and have sync pulses, overwrite tempo drop code times
            %check for actual starting and ending frames of visual stimulus
            %sync pulses are stored in last channel of spike_data
    	    start_eventbin(i) = max(find(all_data.spike_data(num_recorded_spike_channels,:,i) == 1) );
        end
    
        if (stop_code == VSTIM_OFF_CD & size(all_data.spike_data, 1) > 1)
            %if stopcode is stim stop and have sync pulses, overwrite tempo drop code times
            %check for actual starting and ending frames of visual stimulus
            %sync pulses are stored in last channel of spike_data
    	    stop_eventbin(i) = max(find(all_data.spike_data(num_recorded_spike_channels,:,i) == 1) );
        end
    
        if (stop_code == VSTIM_ON_CD & size(all_data.spike_data, 1) > 1)
            %if startcode is stim start and have sync pulses, overwrite tempo drop code times
            %check for actual starting and ending frames of visual stimulus
            %sync pulses are stored in last channel of spike_data
    	    stop_eventbin(i) = min(find(all_data.spike_data(num_recorded_spike_channels ,:,i) == 1) );
        end 
    end

    % here, check for indices lying outside of data range prior to computing spikes
    % can add more conditions
    
    if (start_eventbin(i) + StartOffsetIndex < 0)
        ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
        listtext{length(listtext)+1} = ['ERROR: Start offset smaller than data range for trial ',num2str(i), '.  Truncating start offset.'];
        set(ListHandle, 'String', listtext);
        StartOffsetIndex = 1 - start_eventbin(i);
    end
    if (stop_eventbin(i) + StopOffsetIndex > size(all_data.event_data, 2))
        ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
        listtext{length(listtext)+1} = ['ERROR: Stop offset is greater than data range for trial ',num2str(i), '.  Truncating stop offset.'];
        set(ListHandle, 'String', listtext);
        StopOffsetIndex = size(all_data.event_data, 2) - stop_eventbin(i);
    end   
    if (start_eventbin(i) + StartOffsetIndex > stop_eventbin(i) + StopOffsetIndex)
        ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
        listtext{length(listtext)+1} = ['ERROR: Start offset is set after stop offset.  Stop offset is being reset to match end of data set.'];
        set(ListHandle, 'String', listtext);
        StopOffsetIndex = 1 + size(all_data.event_data, 2) - start_eventbin(i);                    
    end
end

% Now set GUI display to reflect changes in offsets.
StartOffsetHandle = findobj(gcbf, 'Tag', 'StartOffset');			% next four lines added BJP 3/1/00
set(StartOffsetHandle, 'String', num2str(StartOffsetIndex));
StopOffsetHandle = findobj(gcbf, 'Tag', 'StopOffset');
set(StopOffsetHandle, 'String', num2str(StopOffsetIndex));

return;
