function spontRate = FindWinDiscrimSpontRate(handles)
global CHAN2;

% Get the codes times we need to calculate the spontaneous rate.
cTimes = FindSpontaneousCodeTimes(handles);

% Create a vector of all the window discriminator spike times.
eTimes = round([CHAN2(1:length(CHAN2)).time] / handles.chand);

% Go through each set of code times and look for spikes above the threshold.
totalTime = 0;
totalSpikes = 0;
for (i = 1:length(cTimes))    
    spikeIndex = find(eTimes >= cTimes(i).sTime & eTimes <= cTimes(i).eTime);
    totalTime = totalTime + cTimes(i).eTime - cTimes(i).sTime;
    
    % If spikeIndex is empty then go to the next iteration.
    if (isempty(spikeIndex))
        continue;
    end
    
    % Count the number of spikes found in this time segment.
    totalSpikes = totalSpikes +  length(spikeIndex);
end 

totalTime = totalTime / 25000;  % Convert 'totalTime' into seconds.
if totalSpikes == 0 then
    spontRate = 0;
else
    spontRate = round(totalSpikes / totalTime);
end