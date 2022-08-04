%FIND_TIMING_EVENTS Extract event timing information from timing channel,
%assuming one stim-off and one stim-on voltage.
%
% transitions = find_timingchan_events(event_wave, stimpars, sampling_interval_us)
%
% event_wave = signal from stim chan
% stimpars = stimpars struct from stimdata
% sampling_interval_us = sampling_interval in microseconds
%
% %%% Last updated 2018/12/4 on thermaltake
function transitions = find_timingchan_events(event_wave, stimpars, sampling_interval_us)

remove_repeated_values = 0; % toggle.... may need to turn off for parsing resonant galvo scopetiming because events are only 1 interval long

int_wave = round(double(event_wave));  %% round stim wave to nearest integer; use double because uint8 can cause trouble when using diff==0
int_wave = reshape(int_wave, 1, []); % make into row vector
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
    
% remove repeated values
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
if exist('stimpars','var')
    if transitions(2,2)==stimpars.pulse.Vstimoff & transitions(2,1)~=stimpars.pulse.Vstimoff
       warning(['Ignoring stimchan value preceding first stim-off channel value'...
           ' (value = %g from scan 1 to scan %g).'],transitions(2,1),transitions(1,2));
       transitions = transitions(:,2:end); % delete the preceding nonzero column
    end
end
    
% Check that stim events follow the prescribed order of 
% Vstimoff, Vstimon, Vstimoff...Vstimoff
% probably the stimchan value before stimulus presentation began. 
if exist('stimpars','var')
    if rem(length(transitions)-1,2) % if transitions-1 is NOT divisible by 2
        warning('Number of events does not equal 2 * nPresentations + 1')
        input('Press Enter to ignore the last event signal.');
        transitions = transitions(:,1:end-1);
    end
    if any(transitions(2,1:2:end) ~= stimpars.pulse.Vstimoff) ||...    % if Vstimon values are incorrect
           any(transitions(2,2:2:end) ~= stimpars.pulse.Vstimon)      % ... or Vstimoff values are incorrect
       error('Event timing-stimparsing error: events must follow this sequence: [Vstimoff,Vstimon,Vstimoff... VStimoff].')
    end
end

transitions(:,1:2:end) = [];    %% eliminate stim-absent (Vstimoff) events

% Check that response window is shorter than the time between stim onsets. 
% (This step will not the check final stim vs. end of recording.) 
if exist('sampling_interval_us','var') && exist('stimpars','var') && isfield(stimpars,'stimdur')
    stimdur_us = stimpars.stimdur * 1e6; % convert from s to microseconds
    window_scans = stimdur_us / sampling_interval_us; 
    if any(diff(transitions(1,:) < window_scans))
        error('Event timing-stimparsing error: response window to analyze must be shorter than time between stimulus onsets.')
    end
end 
    
transitions = array2table(transitions','VariableNames',{'onset','chan_val','duration'});
