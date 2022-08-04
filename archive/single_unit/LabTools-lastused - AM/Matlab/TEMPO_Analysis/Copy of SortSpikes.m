%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SortSpikes(handles)
%       Sorts ADC Channel spikes based on a specified threshold level and dumps the results in a
%       global structure called 'spsData'.
%
%   @Author: Christopher Broussard
%   @Date:   April 8, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = SortSpikes(handles, isBatch)
global CHAN1 CHAN2 CHAN32 spsData;
index = [];

% Find all the stimulus start codes that resulted in a successful trial.
if (~isBatch)
    set(handles.ProcessText, 'String', 'Finding Good Trials');
else
    disp('Finding Good Trials');
end

for (i = 1:length(CHAN32))
    % If we find a stimuls start code, then we look to see if we find
    % a success code before we get to another start code.
    if(CHAN32(i).mvals == 4)
        j = i + 1;
        while (j <= length(CHAN32) & CHAN32(j).mvals ~= 4)
            if (CHAN32(j).mvals == 12)
                index = [index, i];
                break;
            else
                j = j + 1;
            end
        end
        i = j;
    end
end

% Stuff the spsData struct array with information about each successfull trial.
handles.spikeSet = handles.spikeSet + 1;
spsData(handles.spikeSet).sampleRate = 25000;
spsData(handles.spikeSet).prebuffer = round(handles.PreEventBuffer * 25000);
spsData(handles.spikeSet).postbuffer = round(handles.PostEventBuffer * 25000);
if (get(handles.UpperThreshCheck, 'Value'))
    spsData(handles.spikeSet).upperThreshold = handles.utValue / 5.0 * (2^15 - 1);
else
    spsData(handles.spikeSet).upperThreshold = NaN;
end
if (get(handles.LowerThreshCheck, 'Value'))
    spsData(handles.spikeSet).lowerThreshold = handles.ltValue / 5.0 * (2^15 - 1);
else
    spsData(handles.spikeSet).lowerThreshold = NaN;
end
for (i = 1:length(index))
    if (~isBatch)
        set(handles.ProcessText, 'String', ['Sorting Channel 1 Success Trial ', num2str(i)]);
    else
        disp(['Processing event ', num2str(i)]);
    end
    
    % This is basically descriptive data about the spike area analyzed.
    spsData(handles.spikeSet).spikeInfo(i).startCodeTime = round(CHAN32(index(i)).mark / handles.chand);
    
    spsData(handles.spikeSet).spikeInfo(i).startTime = spsData(handles.spikeSet).spikeInfo(i).startCodeTime - spsData(handles.spikeSet).prebuffer + 1;
    spsData(handles.spikeSet).spikeInfo(i).endTime = spsData(handles.spikeSet).spikeInfo(i).startCodeTime + spsData(handles.spikeSet).postbuffer;
    spsData(handles.spikeSet).spikeInfo(i).wSpikeTimes = [];
    
    % Make sure that the start and end times are within the range of the ADC
    % channel times.
    if (spsData(handles.spikeSet).spikeInfo(i).startTime < 1)
        spsData(handles.spikeSet).spikeInfo(i).startTime = 1;
    end
    if (spsData(handles.spikeSet).spikeInfo(i).endTime > length(CHAN1))
        spsData(handles.spikeSet).spikeInfo(i).endTime = length(CHAN1);
    end
    
    % Store all the event codes for later reference.  Total kludge!!!
    %binCount = (spsData(1).postbuffer + spsData(1).prebuffer) / spsData(1).sampleRate * 1000;
    %slotsperbin = (spsData(1).postbuffer + spsData(1).prebuffer) / binCount;
    %spsData(handles.spikeSet).spikeInfo(i).eventCodes = zeros(1, binCount);
    %mrstart = spsData(handles.spikeSet).spikeInfo(i).startTime * handles.chand;
    %mrend = spsData(handles.spikeSet).spikeInfo(i).endTime * handles.chand;
    %mrsuckass = [CHAN32.mark];
    %mrstupid = find(mrsuckass >= mrstart & mrsuckass <= mrend);
    %a = [CHAN32(mrstupid).mark] / handles.chand;
    %spsData(handles.spikeSet).spikeInfo(i).eventCodes([ceil((a - spsData(handles.spikeSet).spikeInfo(i).startTime + 1) / 25)]) = [CHAN32(mrstupid).mvals];
    
    % Find all the threshold discontinuities in the ADC data.
    upperSpikeIndex = []; lowerSpikeIndex = [];
    if (get(handles.UpperThreshCheck, 'Value'))
        upperSpikeIndex = find(CHAN1(spsData(handles.spikeSet).spikeInfo(i).startTime : spsData(handles.spikeSet).spikeInfo(i).endTime) >= spsData(handles.spikeSet).upperThreshold) + spsData(handles.spikeSet).spikeInfo(i).startTime - 1;
    end
    if (get(handles.LowerThreshCheck, 'Value'))
        lowerSpikeIndex = find(CHAN1(spsData(handles.spikeSet).spikeInfo(i).startTime : spsData(handles.spikeSet).spikeInfo(i).endTime) <= spsData(handles.spikeSet).lowerThreshold) + spsData(handles.spikeSet).spikeInfo(i).startTime - 1;
    end
    
    % Find all the upper threshold spike times.
    if (~isempty(upperSpikeIndex) & get(handles.UpperThreshCheck, 'Value'))
        shiftedData = diff(upperSpikeIndex);
        shiftedData = [2, shiftedData];
        shiftIndex = find(shiftedData > 1);
        
        % Now we actually go through all upper threshold spike areas and find the max.
        for (j = 1:length(shiftIndex))
            % The last element is treated differently
            if (j ~= length(shiftIndex))
                x = [shiftIndex(j) : shiftIndex(j+1) - 1];
                y = upperSpikeIndex(x);
                [n, t] = max(CHAN1(y));
                spsData(handles.spikeSet).spikeInfo(i).pSpikeTimes(j) = upperSpikeIndex(t + shiftIndex(j) - 1);
            else
                x = [shiftIndex(j):length(upperSpikeIndex)];
                y = upperSpikeIndex(x);
                [n, t] = max(CHAN1(y));
                spsData(handles.spikeSet).spikeInfo(i).pSpikeTimes(j) = upperSpikeIndex(t + shiftIndex(j) - 1);
            end
        end
    else
        spsData(handles.spikeSet).spikeInfo(i).pSpikeTimes = [];
    end % End if (~isempty(upperSpikeIndex))
    
    % Find all the lower threshold spike times.
    if (~isempty(lowerSpikeIndex) & get(handles.LowerThreshCheck, 'Value'))
        shiftedData = diff(lowerSpikeIndex);
        shiftedData = [2, shiftedData];
        shiftIndex = find(shiftedData > 1);
        
        % Now we actually go through all upper threshold spike areas and find the max.
        for (j = 1:length(shiftIndex))
            % The last element is treated differently
            if (j ~= length(shiftIndex))
                x = [shiftIndex(j) : shiftIndex(j+1) - 1];
                y = lowerSpikeIndex(x);
                [n, t] = min(CHAN1(y));
                spsData(handles.spikeSet).spikeInfo(i).nSpikeTimes(j) = lowerSpikeIndex(t + shiftIndex(j) - 1);
            else
                x = [shiftIndex(j):length(lowerSpikeIndex)];
                y = lowerSpikeIndex(x);
                [n, t] = min(CHAN1(y));
                spsData(handles.spikeSet).spikeInfo(i).nSpikeTimes(j) = lowerSpikeIndex(t + shiftIndex(j) - 1);
            end
        end
    else
        spsData(handles.spikeSet).spikeInfo(i).nSpikeTimes = [];
    end % End if (~isempty(lowerSpikeIndex))
end % End for (i = 1:length(index))

% Store all the window discriminated spikes only if it was selected to be loaded.
if (handles.maxTime2 > 0)
    eTimes = round([CHAN2(1:length(CHAN2)).time] / handles.chand);
    
    % Set common data about the channel being read in.
    spsData(handles.spikeSet + 1).chanNumber = 2;
    spsData(handles.spikeSet + 1).prebuffer = round(handles.PreEventBuffer * 25000);
    spsData(handles.spikeSet + 1).postbuffer = round(handles.PostEventBuffer * 25000);
    
    for (i = 1:length(index))
        if (~isBatch)
            set(handles.ProcessText, 'String', ['Sorting Channel 2 Success Trial ', num2str(i)]);
        end
        
        % Set descriptive data in spsData.
        spsData(handles.spikeSet + 1).spikeInfo(i).startCodeTime = round(CHAN32(index(i)).mark / handles.chand);
        spsData(handles.spikeSet + 1).spikeInfo(i).startTime = spsData(handles.spikeSet + 1).spikeInfo(i).startCodeTime - spsData(handles.spikeSet + 1).prebuffer + 1;
        spsData(handles.spikeSet + 1).spikeInfo(i).endTime = spsData(handles.spikeSet + 1).spikeInfo(i).startCodeTime + spsData(handles.spikeSet + 1).postbuffer;
        spsData(handles.spikeSet + 1).spikeInfo(i).pSpikeTimes = [];
        spsData(handles.spikeSet + 1).spikeInfo(i).nSpikeTimes = [];
        spsData(handles.spikeSet + 1).spikeInfo(i).wSpikeTimes = eTimes(find(eTimes >= spsData(handles.spikeSet).spikeInfo(i).startTime & eTimes <= spsData(handles.spikeSet).spikeInfo(i).endTime));
    end
end

% Calculate the spontaneous level.
% Go through each set of code times and look for spikes above the threshold.
if (~isBatch)
    set(handles.ProcessText, 'String', 'Calculating Spontaneous Levels');
    cTimes = FindSpontaneousCodeTimes(handles);
    totalTime = 0;
    totalSpikesU = 0; totalSpikesL = 0;
    for (i = 1:length(cTimes))
        % Find upper spontaneous level.
        if (get(handles.UpperThreshCheck, 'Value'))
            spikeIndex = find(double(CHAN1(cTimes(i).sTime : cTimes(i).eTime)) > get(handles.UpperboundSlider, 'Value') /  5.0 * (2^15 - 1));      
            if (~isempty(spikeIndex))
                % Count the number of spikes found in this time segment.
                shiftedData = diff(spikeIndex);
                totalSpikesU = totalSpikesU + length(find(shiftedData > 1)) + 1;
            end
        end
        
        % Find lower spontaneous level.
        if (get(handles.LowerThreshCheck, 'Value'))
            spikeIndex = find(double(CHAN1(cTimes(i).sTime : cTimes(i).eTime)) < get(handles.LowerboundSlider, 'Value') /  5.0 * (2^15 - 1));
            if (~isempty(spikeIndex))
                % Count the number of spikes found in this time segment.
                shiftedData = diff(spikeIndex);
                totalSpikesL = totalSpikesL + length(find(shiftedData > 1)) + 1;
            end
        end
        
        totalTime = totalTime + cTimes(i).eTime - cTimes(i).sTime;
    end  
    
    % Calculate the spike average for both the upper and lower bounds.
    totalTime = totalTime / 25000;  % Convert 'totalTime' into seconds.
    if (get(handles.UpperThreshCheck, 'Value'))
        spikeAverageU = num2str(round(totalSpikesU / totalTime));
    else
        spikeAverageU = 'NULL';
    end
    if (get(handles.LowerThreshCheck, 'Value'))
        spikeAverageL = num2str(round(totalSpikesL / totalTime));
    else
        spikeAverageL = 'NULL';
    end
    
    % Replot data with spike dots.
    handles.plotThreshSpikes = 1;
    [handles.startPoint, handles.endPoint, handles.upperLine, handles.lowerLine] = PlotSpikeData(handles);
    set(handles.upperLine, 'EraseMode', 'xor');
    set(handles.lowerLine, 'EraseMode', 'xor');
    
    set(handles.ProcessText, 'String', 'Sort Complete');
    set(handles.ExportButton, 'Enable', 'on');
    
    % Add this particular spike sort to the Sort List.
    sortList = get(handles.SortList, 'String');
    if (get(handles.UpperThreshCheck, 'Value'))
        ubvalue = num2str(get(handles.UpperboundSlider, 'Value'));
    else
        ubvalue = 'NULL';
    end
    if (get(handles.LowerThreshCheck, 'Value'))
        lbvalue = num2str(get(handles.LowerboundSlider, 'Value'));
    else
        lbvalue = 'NULL';
    end
    s = ['UB: ', ubvalue, ', US: ', spikeAverageU, ' -- LB: ', lbvalue, ', LS: ', spikeAverageL];
    sortList{handles.spikeSet} = s;
    set(handles.SortList, 'String', sortList);
    set(handles.SortList, 'Value', handles.spikeSet);
    
    % Enable the remove and clear buttons.
    set(handles.RemoveSortButton, 'Enable', 'on');
    set(handles.ClearSortListButton, 'Enable', 'on');
end % End of if (~isBatch)

return;
