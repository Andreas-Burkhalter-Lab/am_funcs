%Editted by RLS on 11-25-07 to average over the unstimulated time only
%Editted by RLS on 12-8-07 from calculations in terms of seconds to
%calculations in terms of milliseconds
function [response_time,EOB,P]=burst(InTrain, StartT, StopT)
%NOTE: Calls POISSCDF function in the STATISTICS toolbox
%Adapted from Hanes et al 1995 method for determining when a cell becomes
%active
%
%Input Arguments:
%'InTrain' is the timestamps of all the spikes recorded for a given cell,
%typically for only the trials with the cell's preferred stimulus value
%'StartT' is the timestamp for when the stimulus was actually turned on
%'StopT' is not necessarily when the stimulus presentation ended, but rather
%when the recording for the given trial ended.
%
%Returns:
%'response_time' is how long it takes the cell to become active after the
%stimulus is presented
%'EOB' is the time from stimulus presentation to the end of the cell's
%active burst; this value is not very important, but it can be insightful if
%the function isn't working correctly
%P is the probability that the spike events that fell between
%'response_time' and 'EOB' occurred due to chance assuming that the 
%interspike intervals follow a poisson distribution
%
%NOTE: response_time value of 9999 indicates that the cell did not become
%active during stimulus presentation.  Also, be careful using this with
%data that was recorded before the sync pulses were used.
%
%Created by Chris Andrews 8/13/07

%Define the minimum number of spikes in a burst, so that a very small
%number of closely spaced spikes (2 or 3) occurring at an arbitrary time
%are not misinterpreted as a burst
min_spikes_in_burst = 5;

%Sort spikes from first to last in case they are not already sorted
InTrain = sort(InTrain);

%Define the time that the stimulation starts after the recording begins
%StimStartT = 500;

%Define the time that the stimulation stops during the recording
%StimStopT = 2500;

%While scrolling through all the spike timestamps, we will count how many
%spikes occur during the times before and after stimulation
%implementation
m = 1;
unstim_spikes = 0;
while m <= length(InTrain);
    if InTrain(m) <= 500
        unstim_spikes = unstim_spikes + 1;
        m = m + 1;
    elseif InTrain(m) >= 2500;
        unstim_spikes = unstim_spikes + 1;
        m = m + 1;
    else
        m = m + 1;
    end
end
    
%The mean firing rate must be determined from the entire trial, not just
%the active part, so that this statistical method will work
%MeanRate = length(InTrain)/((StopT-1)*0.001);%spikes per second

%The mean firing rate is now determined only from the time before and after
%the stimulation is being implemented.
MeanRate = unstim_spikes/1000;%spikes/ millisecond (.5 seconds before + .5 seconds after)

%Find the first two spikes whose rate exceeds the null rate
start_flag = 0;
n = 1;
n_spikes = length(InTrain);
while start_flag == 0;
    spike1 = InTrain(n);
    spike2 = InTrain(n+1);
    %Calculate time between spikes in milliseconds
    delta_t = (spike2 - spike1);
    rate = 2/(delta_t);
    if rate > MeanRate
        start_flag = 1;
        start_of_burst = n;
        first_index = n + 1;
    elseif n >= (n_spikes - 2);
        disp('Error: No two spikes exceeded the null rate')
        break
    else
        n = n + 1;
    end
end

%Using the two spikes found above as a starting point, and ensuring that
%the minimum number of spikes in a burst is satisfied, calculate the
%probability that the spikes between start_of_burst and index occurred by
%chance.  Advance all the way to the last spike in the spike train. The end
%of the burst will be the spike which gives the lowest probability of
%chance occurrence.
P_threshold = 0.05;
min_P = P_threshold;
while (exist('end_of_burst') == 0);

    for index = first_index+min_spikes_in_burst-2:n_spikes;
        delta_t = (InTrain(index) - InTrain(start_of_burst));
        spikes_in_burst = (index - start_of_burst + 1);
        lambda = MeanRate*delta_t;
        P = 1 - poisscdf(spikes_in_burst,lambda);
        if (P < min(P_threshold,min_P) && InTrain(index) < StopT)
            end_of_burst = index;
            min_P = P;
            break
        end
    end

    if (exist('end_of_burst') == 0)
        first_index = first_index+1;
    end
    
    if (first_index == length(InTrain))
        end_of_burst = length(InTrain)
        break
    end

end

%Now that the end of the burst is determined, advanced backwards in time to
%find the spike which gives the lowest probability that the burst occurred
%by chance.  This spike corresponds to the time that the cell became
%active.
min_P = P_threshold;
for index = end_of_burst-min_spikes_in_burst:-1:start_of_burst;
    delta_t = (InTrain(end_of_burst) - InTrain(index));
    spikes_in_burst = (end_of_burst - index + 1);
    lambda = MeanRate*delta_t;
    P = 1 - poisscdf(spikes_in_burst,lambda);
    if (P < min(P_threshold,min_P) && InTrain(index) > StartT)
        onset_time = index;
        min_P = P;
    end
end

%If the cell did not become active, set onset_time to 9999 so that it can
%be identified later.
if (exist('onset_time') == 0)
    onset_time = 9999;
end

if onset_time~=9999;
    response_time = (InTrain(onset_time)-StartT)/1000;
    EOB = (InTrain(end_of_burst)-StartT)/1000;
    P = min_P;
else
    response_time = 9999;
    EOB = 9999;
    P = 9999;
end