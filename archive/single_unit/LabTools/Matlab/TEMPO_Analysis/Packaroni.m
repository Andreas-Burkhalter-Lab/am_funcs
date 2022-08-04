%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function [good_data] = Packaroni(good_data, fileName, wSpikesOnly, shiftValue)
%       Packaroni takes the good_data data structure an stuffs it with the spiketimes5
%
%   @Author: Christopher Broussard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [good_data, outputChan] = Packaroni(good_data, fName, wSpikesOnly, shiftValue)

TEMPO_Defs;

% Load in the spike data from spikesort2.
load(fName);

global DescriptionsExist;
DescriptionsExist = 0;

% If spsData is empty just return.
if (isempty(spsData))
    return;
end

% Find the next channel # to store data in.
chanNumber = size(good_data.spike_data, 1);

% The number of bins and slots per bin
binCount = (spsData(1).postbuffer + spsData(1).prebuffer) / spsData(1).sampleRate * (good_data.htb_header{SPIKE_DB}.speed_units/good_data.htb_header{SPIKE_DB}.speed);
slotsperbin = (spsData(1).postbuffer + spsData(1).prebuffer) / binCount;
n = zeros(1, binCount);

% If all we want is the window discriminator spikes, then don't pack the sorted
% spikes.
if (wSpikesOnly == 1)
    index = length(spsData);
    for (j = 1:length(spsData(index).spikeInfo))
        % Bin all the spike times.
        edges = [spsData(index).spikeInfo(j).startCodeTime - spsData(index).prebuffer + 1 : slotsperbin : spsData(index).spikeInfo(j).startCodeTime + spsData(index).postbuffer + slotsperbin];
        combinedSpikes = sort([spsData(index).spikeInfo(j).pSpikeTimes, spsData(index).spikeInfo(j).nSpikeTimes, spsData(index).spikeInfo(j).wSpikeTimes]);
        if (isempty(combinedSpikes))
            n = zeros(1, binCount);   
        else
            n = histc(combinedSpikes, edges);
            n = n(1:length(n) - 1);
        end
        
        % Put the bins in good_data.             
        good_data.spike_data(chanNumber + 1, :, j) = n;
        outputChan = chanNumber + 1;
    end    
    % Store sorted spike channel tag.
    if (exist('descExists') > 0)
        DescriptionsExist = 1;
        good_data.desc{chanNumber + 1} = spsData(index).desc;
    else
        good_data.desc{chanNumber + 1} = [];
    end
else  
    % Go through all the data elements of spsData and stuff the spike times
    % into good_data.
    for (i = 1:length(spsData) - wdata) %loop per spike template
        for (j = 1:length(spsData(i).spikeInfo)) %loop per trial
            % Bin all the spike times.
            edges = [spsData(i).spikeInfo(j).startCodeTime - spsData(i).prebuffer + 1 : slotsperbin : spsData(i).spikeInfo(j).startCodeTime + spsData(i).postbuffer + slotsperbin];
            combinedSpikes = sort([spsData(i).spikeInfo(j).pSpikeTimes, spsData(i).spikeInfo(j).nSpikeTimes, spsData(i).spikeInfo(j).wSpikeTimes]);
            if (isempty(combinedSpikes))
                n = zeros(1, binCount);   
            else
                n = histc(combinedSpikes, edges);
                n = n(1:length(n) - 1);
            end
            
            % Shift the data according to the shiftValue.
            if (shiftValue > 0)
                n = [zeros(1, shiftValue), n];
                n = n(1:length(n) - shiftValue);
            elseif (shiftValue < 0)
                n = [n, zeros(1, abs(shiftValue))];
                n = n(abs(shiftValue) + 1 : length(n));
            end
            
            % Put the bins in good_data.
            good_data.spike_data(chanNumber + i, :, j) = n;
            outputChan(i) = chanNumber + i;
        end
        
        % Store sorted spike channel tag.
        if (exist('descExists') > 0)
            DescriptionsExist = 1;
            good_data.desc{chanNumber + i} = spsData(i).desc;
        else
            good_data.desc{chanNumber + i} = [];
        end
    end
end % End if (wSpikesOnly == 1).

return;
