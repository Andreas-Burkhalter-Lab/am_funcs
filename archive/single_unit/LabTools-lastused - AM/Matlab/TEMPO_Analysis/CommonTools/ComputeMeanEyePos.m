%----------------------------------------------------------------------------------------------------------
%-- ComputeMeanEyePos.m: This function computes the average position for each eye channel and each
%--	trial, over a period of time specified by start_code and stop_code, which may occur at different
%--	times on different trials.  JDN, 3/6/00
%----------------------------------------------------------------------------------------------------------
function eye_positions = ComputeMeanEyePos(data, n_trials, start_code, stop_code, StartOffsetBin, StopOffsetBin);

TEMPO_Defs;	%some defines that we'll need

DEBUG = 0;

%first, we'll need the bin_width (in sec) of our eye channels + events log
h = data.htb_header{EYE_DB};	%for convenience
eye_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

% compute indexing in eye database, start index at 1
StartOffset = 1 + floor (StartOffsetBin*(event_bin_width/eye_bin_width));
StopOffset = 1 + floor (StopOffsetBin*(event_bin_width/eye_bin_width));


%Here, I count eye positions over the period from start_index to stop_index for every trial
%In many cases, start_index and stop_index would be the same for each trial, and the loops would not be needed
%But, to be general, I am doing it this way so that each trial could have different start- and stop_ indices
for j = 1:data.htb_header{EYE_DB}.nchannels		%for each eye channel (=4)
    if (DEBUG)
        figure;
        hold on;
    end
    for i = 1:n_trials		%for each trial
        start_eventbin = find(data.event_data(1,:,i) == start_code);   
        stop_eventbin = find(data.event_data(1,:,i) == stop_code);   
        
        if ( isempty(start_eventbin) | isempty(stop_eventbin) )
            start_eventbin = NaN; 
            stop_eventbin = NaN;
        end
        
        % convert to eye bins
        start_eyebin = NaN; stop_eyebin = NaN;
        start_eyebin = floor (start_eventbin*(event_bin_width/eye_bin_width)) + StartOffset;
        stop_eyebin = floor (stop_eventbin*(event_bin_width/eye_bin_width)) + StopOffset;
        
        if ( isnan(start_eyebin) | isnan(stop_eyebin) )
            eye_positions(j,i) = NaN; 
        else
            %sum the eye positions for this trial
            count = sum(data.eye_data(j,start_eyebin:stop_eyebin,i), 2);
            %divide the summed eye positions by the # of bins to get average position
            eye_positions(j,i) = count / (stop_eyebin - start_eyebin); 
            
            if (DEBUG)
                plot(data.eye_data(j,start_eyebin:stop_eyebin,i), 'k-');
            end
        end
    end
    if (DEBUG)
        hold off;
    end
end

return;
