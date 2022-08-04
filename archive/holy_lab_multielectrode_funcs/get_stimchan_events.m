%GET_STIMCHAN_EVENTS Extract event timing information from stim channel,
%assuming one stim-off and one stim-on voltage.
%%% Last updated 1/7/15 on msi

%% added removal of repeated values.... came up in 10/2/15 recording
remove_repeated_values = 1; % toggle

%% add integration of mouse_movements here

if ~exist('stimrec','var')
    if exist('stimrec_temp','var') && ~isempty('stimrec_temp')
        fprintf(['Variable ''stimrec'' not found in stimulus log file (stim presentation ',...
            'probably aborted before completion). \n     Using ''stimrec_temp'' instead.\n'])
        stimrec = stimrec_temp; % copy-on-write
    else
        error('Did not find ''stimrec'' or ''stimrec_temp'' in stimulus log file.')
    end
end

int_wave = round(merec_obj(stim_chan, 1:merec_obj.nscans));  %% round stim wave to nearest integer
transitions(1,:) = [1, find(~diff(int_wave)==0) + 1]; % event times, starting with initial state (1)... +1 is to find post-transition time
transitions(2,:) = int_wave(transitions(1,:)); %% new stim channel values after transitioning from one value to another   
transitions(3,:) = [diff(transitions(1,:)) NaN]; % event durations  (ignore last nonstim event duration)

% Delete all event transitions immediately following the preceding transition, 
% in case it takes more than one timestep to complete a transition. 
trans_diff = [diff(transitions(1,:)) 0]; %% get time length in time steps between events, append 0 for last timestep
countback = 0;
for i = 1:size(transitions,2)
    if trans_diff(i)==1 %% if this transition is 1 timestep before the next transition
        countback = countback+1; 
    elseif countback  %% if this transition is 1 step away from the preceding but not the following transition
        transitions(1,i-countback+1:i) = -1;        %% mark problematic transitions with -1 for deletion
        transitions(2,i-countback) = transitions(2,i); % use the final stimchan value in this set of adjacent transitions
        transitions(3,i-countback) = transitions(3,i)+countback; %% correct the duration of nondeleted event
        countback = 0;
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

% If the initial pre-stim Vstimoff value is preceded by another stimchan
% value, give a warning but delete the first value. The first value was
% probably the stimchan value before stimulus presentation began. 
if transitions(2,2)==pulsepars.Vstimoff & transitions(2,1)~=pulsepars.Vstimoff
   warning(['Ignoring stimchan value preceding first stim-off channel value'...
       ' (value = %g from scan 1 to scan %g).'],transitions(end,1),transitions(1,2));
   transitions = transitions(:,2:end); % delete the preceding nonzero column
end

% Check that stim events follow the prescribed order of 
% Vstimoff, Vstimon, Vstimoff...Vstimoff
if rem(length(transitions)-1,2) % if transitions-1 is NOT divisible by 2
    warning('Number of events does not equal 2 * nPresentations + 1')
    input('Press Enter to ignore the last event signal.');
    transitions = transitions(:,1:end-1);
end
if any(transitions(2,1:2:end) ~= pulsepars.Vstimoff) ||...    % if Vstimon values are incorrect
       any(transitions(2,2:2:end) ~= pulsepars.Vstimon)      % ... or Vstimoff values are incorrect
   error('Event timing-parsing error: events must follow this sequence: [Vstimoff,Vstimon,Vstimoff... VStimoff].')
end

transitions(:,1:2:end) = [];    %% eliminate stim-absent (Vstimoff) events

% Check that response window is shorter than the time between stim onsets. 
% (This will not the check final stim vs. end of recording.) 
if any(diff(transitions(1,:) < window_scans))
    error('Event timing-parsing error: response window to analyze must be shorter than time between stimulus onsets.')
end

% Check for stim file/.merec file trial number mismatch.
%%%%%%%% maybe we should list the actual number of trials acquired for a
%%%%%%%% stim param combination, because if trials are missing then trials
%%%%%%%% per combination are probably not equal
if size(transitions,2) < size(stimrec,1) % if we have electrode data for fewer trials than we have stim data for
    warning(['Only %g trials were found in electrode data out of %g trials ',...
        'listed in the stimulus file. Trials without electrode data will be ignored.'],...
        size(transitions,2), size(stimrec,1));
    stimrec = stimrec(1:size(transitions,2),:); % remove trial descriptions for which we don't have electrode data
elseif size(transitions,2) > size(stimrec,1) % if we are missing stim data for some trials
        warning(['Only %g trials were listed in stimulus file data out of %g trials',...
        'found in the electrode recording file. Trials without stimulus parameter data will be ignored.'],...
        par_sets, size(transitions,2)); % because trialdata comes from stimrec, extra trials will be passively ignored
end

trialdata = stimrec; % rename with copy-on-write
clear stimrec; 
% trialdata.diam_pix = []; % remove unneeded fields 
% trialdata.sfreq_pix = []; % remove unneeded fields  
trialdata.onset = (transitions(1,:))'; % stim onset time in scans
trialdata.dur = (transitions(3,:))'; % duration in scans
clear transitions trans_diff count
save([getfname(resp_file) '.trialdata'],'trialdata'); % save in case we want to run analyze_ss_responses again