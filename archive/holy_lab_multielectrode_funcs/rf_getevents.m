%RF_GETEVENTS Get event timing info from the stim channel for
%analyze_rf_responses.
%%% Last updated 1/30/16 on vivid

%% added removal of repeated values.... came up in 10/2/15 recording with preamp
remove_repeated_values = 1; % toggle

%% Extract and check stim channel signals
% This script assumes the following: V is held at pulse.Vbase when no stimuli
% are present, then shifts to a negative integer at stimulus onset, then
% shifts to a positive integer and is held there until stim offset, at which
% time V returns to Vbase; also, that no stim are present at t=0 and Vbase
% is an integer. 
% 'stim_info' contains event onsets in row 1, event
%   durations in row 2, and compressed event labels in row 3. 

stimwave = merec_obj([stim_chan], [1:merec_obj.nscans]);
window_scans = resp_window * merec_obj.scanrate;    % convert window size from seconds to scans

int_wave = round(stimwave);  %% round stim wave to nearest integer
transitions(1,:) = [1, find(~diff(int_wave)==0) + 1]; % event timings, starting with initial state (1)... +1 is to find post-transition time
transitions(2,:) = int_wave(transitions(1,:)); %% get stim channel values after transitioning from one value to another

% Delete all event transitions immediately following the preceding transition, 
% in case it takes more than one timestep to complete a transition. 
trans_diff = [diff(transitions(1,:)) 0]; %% get time length in time steps between events, append 0 for last timestep
count = 0;
for i = 1:size(transitions,2)
    if trans_diff(i)==1 %% if this transition is 1 timestep before the next transition
        count = count+1; 
    elseif count  %% if this transition is 1 step away from the preceding but not the following transition
        transitions(1,i-count+1:i) = -1;        %% mark problematic transitions with -1 for deletion
        transitions(2,i-count) = transitions(2,i); % use the final stimchan value in this set of adjacent transitions
        count = 0;
    end
end
transitions = transitions(:,transitions(1,:)>0);  %% eliminate all transitions one timestep after the preceding transition

%% remove repeated values.... why are these appearing?
if remove_repeated_values
   countback = 0;
   total_countback = 0;
   for i = 1:length(transitions)
        if i == length(transitions) % if last event
            if countback > 0 % if this is the end of a repeated sequence
               transitions(2,i-countback+1 : i) = inf; % mark all in the repeated sequence except the first for deletion  
            end
        elseif transitions(2,i+1) == transitions(2,i) % if the next value is the same as this value
            countback = countback+1;
            total_countback = total_countback+1; 
        elseif countback > 0 % if this is the end of a repeated sequence
            transitions(2,i-countback+1 : i) = inf; % mark all in the repeated sequence except the first for deletion     
            countback = 0; % reset
        end
   end
   trans_inf = transitions(2,:) == inf;
   transitions(:,trans_inf) = [];
   transitions(3,:) = [diff(transitions(1,:)) NaN]; % correct event durations  (ignore last nonstim event duration)
   if total_countback > 0
       fprintf('Removed %g contiguous repeated event values; %g events remain.\n',total_countback,size(transitions,2));
   end
end

% If the initial pre-stim Vbase value is preceded by another stimchan
% value, give a warning but delete the first value. The first value was
% probably the stimchan value before stimulus presentation began. 
if transitions(2,2)==pulse.Vbase & transitions(2,1)~=pulse.Vbase % if 2nd is vbase and 1st is not vbase
   warning(['Ignoring stimchan value preceding first stim-off channel value'...
       ' (value = %g from scan 1 to scan %g).'],transitions(end,1),transitions(1,2));
   transitions = transitions(:,2:end); % delete the preceding nonzero column
end

% Check that stim events follow the prescribed order of of Vbase,negative,positive,Vbase....positive,Vbase
if rem(length(transitions)-1,3) % if 'transitions'-1 is NOT divisible by 3
    error('Event timing-parsing error: number of events must equal 3 * nUniqueStim + 1')
elseif any(transitions(2,1:3:end) ~= pulse.Vbase)pulse ||...    % if Vbase values are inccorrect
       any(transitions(2,2:3:end) >= 0) || ...       % ... or negative values are incorrect
       any(transitions(2,3:3:end) <= 0)              % .... or positive values are incorrect
   error('Event timing-parsing error: events must follow this sequence: [Vbase, -, +, Vbase, -.... +, Vbase].')
end

% Compress labels: reverse sign of negative events and multiply by max_vsteps, then add the following positive event, 
% then delete the positive events (positive events only encode information about the preceding negative event). 
transitions(2,2:3:end) = -max_vsteps*transitions(2,2:3:end) + transitions(2,3:3:end); % negative event = tens place, positive event = ones place
transitions(:,3:3:end) = []; % eliminate positive-event columns




%% would probably be faster and simpler to get rid of the 'event_bins' variable and just use 'trialdata_rf' instead


% Convert to event-block form as function output: first row is stim onset, 
% second row is duration, third row is event label. 
event_bins = [transitions(1,:); [diff(transitions(1,:)) 0]; transitions(2,:)];%add 2nd row: event durations (ignore last nonstim event duration)
event_bins(:,1:2:end) = [];    %% eliminate stim-absent (Vbase) events

% Check that response window is shorter than the time between stim onsets. (This will not the check final stim vs. end of recording.) 
if any(diff(event_bins(1,:) < window_scans))
    error('Event timing-parsing error: response window to analyze must be shorter than time between stimulus onsets.')
end


% Checks: is #events correct, is each row/column combo found once per
% iteration, and are iterations identical? 
%%%%%%%%%%%%%%%% maybe add a check that stim duration is the expected length (during event parsing section)    
nevents = size(event_bins,2);
if stimpars_rf.nAngles*(stimpars_rf.rows * stimpars_rf.columns) ~= nevents % check that event_bins is the expected size
    error(['Event timing-parsing error: number of detected stim events (%g)',...
        ' does not equal of iterations (%g) x specified rows (%g) x specified columns (%g) = %g.'],...
        size(event_bins,2), stimpars_rf.nAngles, stimpars_rf.rows, stimpars_rf.columns,...
        stimpars_rf.nAngles*stimpars_rf.rows*stimpars_rf.columns);
else  % if the number of events is correct
    first_iter = event_bins(:, 1 : (size(event_bins,2)/stimpars_rf.nAngles)); % first iteration of stimuli
    for row = stimpars_rf.rows
        for column = stimpars_rf.columns
            if numel(find( 10*row + column == first_iter(3,:))) ~= 1 % if this row/column combination isn't represented exactly once
                error(['Event timing-parsing error: in each iteration, each row/column combination',...
                    'must be represented exactly once on the stim channel.'])
            end
        end
    end
    if event_bins(3,:) ~= repmat(first_iter(3,:), 1, stimpars_rf.nAngles); % if each iteration isn't identical
        error('Event timing-parsing error: each iteration must be identical')
    end
end

% Put event_bins into a dataset to make easier to read and save event data.
% save in case we want to run rf_main again
events_per_iter = nevents/stimpars_rf.nAngles;
iterLabels = [kron(1:stimpars_rf.nAngles,ones(1,events_per_iter))]';
trialdata_rf = dataset(event_bins(1,:)',event_bins(2,:)',event_bins(3,:)',iterLabels,... 
    'VarNames',{'onset','dur','rowcolumn','iteration'});
save([getfname(resp_file_rf) '.trialdata'],'trialdata_rf','event_bins','window_scans','nevents','first_iter'); 



