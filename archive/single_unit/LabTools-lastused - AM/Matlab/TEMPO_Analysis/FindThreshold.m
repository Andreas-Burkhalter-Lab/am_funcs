% --------------------------------------------------------------------------
%   This is a recursive function that tries to find a threshold value
%   that closely matches a specified spontaneity value.  The polarity
%   parameter lets you choose whether you're looking at positive or
%   negative numbers.
% --------------------------------------------------------------------------
function threshValue = FindThreshold(cTimes, rangeDivide, handles, tValue, spontValue, polarity, isBatch)
    global CHAN1;
    totalTime = 0;
    totalSpikes = 0;
    
    if (~isBatch)
        set(handles.ProcessText, 'String', ['Checking threshold: ', num2str(tValue)]);
    end
    
    % Go through each set of code times and look for spikes above the threshold.
    for (i = 1:length(cTimes))
        spikeIndex = find(double(CHAN1(cTimes(i).sTime : cTimes(i).eTime)) * polarity >= (tValue / 5.0 * (2^15 - 1)));
        totalTime = totalTime + cTimes(i).eTime - cTimes(i).sTime;
        
        % If spikeIndex is empty then go to the next iteration.
        if (isempty(spikeIndex))
            continue;
        end
        shiftedData = diff(spikeIndex);
        
        % Count the number of spikes found in this time segment.
        totalSpikes = totalSpikes + length(find(shiftedData > 1)) + 1;
    end  
    
    totalTime = totalTime / 25000;  % Convert 'totalTime' into seconds.
    spikeAverage = round(totalSpikes / totalTime);
    
    % If the spikeAverage is equal to the spontValue, then we can return the ideal
    % threshold.  Otherwise, keep searching.
    % If our spike average is < spontValue, then we need to increase our tValue.
    % If it's >, then we need to decrease tValue.
    if (spikeAverage == spontValue)
        threshValue = tValue;
    elseif (spikeAverage < spontValue)
        threshValue = FindThreshold(cTimes, rangeDivide + 1, handles, tValue - (handles.rangeMax / 2^rangeDivide), spontValue, polarity, isBatch);
    else
        threshValue = FindThreshold(cTimes, rangeDivide + 1, handles, tValue + (handles.rangeMax / 2^rangeDivide), spontValue, polarity, isBatch);
    end
    
    return;
