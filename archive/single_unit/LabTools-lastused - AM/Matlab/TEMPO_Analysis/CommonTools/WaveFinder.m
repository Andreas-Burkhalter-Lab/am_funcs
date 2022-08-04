function waveData = WaveFinder(fileName, htbFileName)

% Input error checking.
if (exist('fileName') == 0)
    error('Please specify an SMR file.');
end
if (exist('htbFileName') == 0)
    error('Please specify an htb file.');
end

fprintf('Starting WaveFinder\n');

% Load Tempo info
[pretime, posttime, config] = LoadTempoFileData(htbFileName);

% Open the SMR file.
handle = cedfunction('SonOpenOldFile', fileName, 1);
if (handle < 0)
    return;
end

infoStruct.smrFile = fileName;
infoStruct.tFile = htbFileName;
infoStruct.preTime = pretime;
infoStruct.postTime = posttime;
infoStruct.handle = handle;
infoStruct.chand = cedfunction('SonChanDivide', handle, 0);
infoStruct.config = config;
% Get the start code times.
fprintf('Finding start code times.\n');
startCodeTimes = FindStartCodeTimes(infoStruct);
infoStruct.startCodeTimes = startCodeTimes;

% Find all ADCMark channels.
waveChans = [];
for (i = 0:31)
    if (cedFunction('SonChanKind', handle, i) == 6)
        waveChans = [waveChans, i];
    end
end



% Grab the ADCMark code and times.
waveData = [];
index = 1;
for (i = waveChans)
    fprintf('Processing wavemark channel: %d\n', i + 1);  %CED channels are 0 - 31, screen 1-32
    
    % Get the channel maxtime.
    maxtime = cedFunction('SonChanMaxTime', handle, i);
    
    % Get the marker data.
    [num, x] = cedFunction('SonGetMarkData', handle, i, 100000, 0, maxtime);
    
    % Find all unique marker values in the channel.
    uniqueVals = munique([x.mvals]')';
    
    % noise interpreted as marker -18, unclassified = 0.  Do not process
    % these - BJP 10/21/03
    
    uniqueVals = uniqueVals(uniqueVals > 0);
    for (j = uniqueVals)
        fprintf('Processing marker code: %d\n', j);
        
        waveData(index).sampleRate = 25000;
        waveData(index).prebuffer = pretime;
        waveData(index).postbuffer = posttime;
        waveData(index).desc = ['Sorted Channel: ', num2str(i + 1), ', WaveMark: ', num2str(j)];
        
        markTimes = [x(find([x.mvals] == j)).mark] / infoStruct.chand; %in CED units
        
        % Separate the marker times into trials.
        for (k = 1:length(startCodeTimes))
            waveData(index).spikeInfo(k).startCodeTime = startCodeTimes(k);
            waveData(index).spikeInfo(k).startTime = startCodeTimes(k) - waveData(index).prebuffer + 1;
            waveData(index).spikeInfo(k).endTime = startCodeTimes(k) + waveData(index).postbuffer;
            waveData(index).spikeInfo(k).nSpikeTimes = [];
            waveData(index).spikeInfo(k).wSpikeTimes = [];
            
            % Find spikes within the trial.
            mindex = find((markTimes >= waveData(index).spikeInfo(k).startTime) & (markTimes <= waveData(index).spikeInfo(k).endTime));
            waveData(index).spikeInfo(k).pSpikeTimes = markTimes(mindex);
        end
        
        index = index + 1;
    end
end

%create the path and filename of output file
dataPath = [infoStruct.smrFile(1:length(infoStruct.smrFile))];
i = size(dataPath,2) - 1;
while (dataPath(i) ~= '\')	%Analysis directory is one branch below Raw Data Dir
     i = i - 1;
end   
dataPath = dataPath(1:i);
i = size(dataPath,2) - 1;
while (dataPath(i) ~= '\')	%Analysis directory is one branch below Raw Data Dir
     i = i - 1;
end   
dataPath = [dataPath(1:i) 'SortedSpikes\'];
        
dataFile = htbFileName;
i = size(dataFile,2) - 1;
while dataFile(i) ~= '.'
    i = i - 1;
end
j = size(dataFile,2) - 1;
while dataFile(j) ~= '\'	%Analysis directory is one branch below Raw Data Dir
    j = j - 1;
end   
dataFile = dataFile(j + 1:i - 1);

OutputDataFile = [dataPath dataFile '.mat'];

% Check to see if it exists.
saveOpts = '';
descExists = 1;
if (exist(OutputDataFile, 'file') == 0)
    fprintf('Writing new data file.\n');
    spsData = waveData;
    
    [winflag, discrimData] = LoadWindowDiscrimData(infoStruct);
    if (winflag)
        spsData(length(spsData)+1) = discrimData;
        wdata = 1;
    else
        wdata = 0;
    end
    
    saveOpts = ' wdata descExists';
else
    fprintf('Adding to existing data file.\n');
    
    % Load 'spsData' and append the new data.
    load(OutputDataFile);
    eval('offset = length(spsData);');
    saveOpts = ' descExists -append';
    
    % Store the window discriminator data.
    discrimData = spsData(offset);
    
    % Append the new data, overwriting the old window discriminator
    % channel.
    for (i = 1:length(waveData))
        spsData(i+offset-1) = waveData(i);
    end
    
    % Re-add the window discriminator channel at the end.
    spsData(length(spsData)+1) = discrimData;
end

eval(['save ', OutputDataFile,  ' spsData', saveOpts, ';']);

fprintf('Finished\n');

cedfunction('SonCloseFile', handle);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads window discriminator data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wnumpts, winData] = LoadWindowDiscrimData(infoStruct)
% Load ADC Channel 2 into memory (Window Discriminator).
maxTime2 = cedfunction('SonChanMaxTime', infoStruct.handle, 1);
if (maxTime2 > 0)
    [num, CHAN2] = cedfunction('SonGetEventData', infoStruct.handle, 1, 1000000, 0, maxTime2);
    
    winData.sampleRate = 25000;
    winData.prebuffer = infoStruct.preTime;
    winData.postbuffer = infoStruct.postTime;
    winData.desc = ['Windows Discriminator Channel'];
    
    markTimes = round([CHAN2(1:length(CHAN2)).time] / infoStruct.chand); %in CED units
    
    % Separate the marker times into trials.
    for (k = 1:length(infoStruct.startCodeTimes))
        winData.spikeInfo(k).startCodeTime = infoStruct.startCodeTimes(k);
        winData.spikeInfo(k).startTime = infoStruct.startCodeTimes(k) - winData.prebuffer + 1;
        winData.spikeInfo(k).endTime = infoStruct.startCodeTimes(k) + winData.postbuffer;
        winData.spikeInfo(k).nSpikeTimes = [];
        winData.spikeInfo(k).pSpikeTimes = [];
        
        % Find spikes within the trial.
        mindex = find((markTimes >= winData.spikeInfo(k).startTime) & (markTimes <= winData.spikeInfo(k).endTime));
        winData.spikeInfo(k).wSpikeTimes = markTimes(mindex);
    end
    
    wnumpts = 1;
else
    wnumpts = 0;
    winData = [];
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns array of start code times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startCodeTimes = FindStartCodeTimes(infoStruct)

startCodeTimes = [];

% Load the marker channel.
maxTime32 = cedfunction('SonChanMaxTime', infoStruct.handle, 31);
[num, chan32] = cedfunction('SonGetMarkData', infoStruct.handle, 31, 1000000, 0, maxTime32);

%klug because my marker codes seem to be close to doubled _BJP
% if ~strcmp(infoStruct.config, 'DotsDisc.pcf' )
%     maxTime32 = maxTime32 / 2;
%     for i = 1: length(chan32)
%         chan32(1).mark = chan32(1).mark / 2;
%     end
% end     
    
% Fix the marker channel.
[fixable, errsFound, chan32] = FixData(chan32, infoStruct);
if (fixable == 0)
    cedfunction('SonCloseFile', infoStruct.handle);
    error('Could not fix SMR file, closing file!');
end

% Find start codes.
index = [];
for (i = 1:length(chan32))
    % If we find a stimuls start code, then we look to see if we find
    % a success code before we get to another start code.
    if(chan32(i).mvals == 4)
        j = i + 1;
        while (j <= length(chan32) & chan32(j).mvals ~= 4)
            if (chan32(j).mvals == 12)
                index = [index, i];
                break;
            else
                j = j + 1;
            end
        end
        i = j;
    end
end

for (i = index)
    startCodeTimes = [startCodeTimes, round(chan32(i).mark / infoStruct.chand)];
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Determines the pre and post event buffer time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [preEventBuffer, postEventBuffer, config] = LoadTempoFileData(tfile)

fid = htbOpen(tfile);     % Open the htb file.
ndbs = htbCount(fid);     % Find out the # of databases it holds.

% Find the Events database and get the times.
for (i = 1:ndbs)
    hd = htbGetHd(fid, i);   % Get the database header.
    
    if (strcmp(hd.title, 'Events'))
        hertz = hd.speed_units / hd.speed;    % see p366 in TEMPO v9 manual
        binwidth = (hd.skip + 1) / hertz;
        epoch_start = hd.offset * binwidth;
        epoch_period = hd.period * binwidth;
        preEventBuffer = epoch_start;
        postEventBuffer = -epoch_start + epoch_period;
        break;
    end
end

% Convert to CED units.
preEventBuffer = round(preEventBuffer * 25000);
postEventBuffer = round(postEventBuffer * 25000);
config = hd.cfg_file;
htbClose(fid);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Finds and fixes problems with SMR data that disagrees with Tempo
%       data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fixable, isRecovered, chan32] = FixData(chan32, infoStruct)
fixable = 1;
isRecovered = 0;
index = [];

% Check to see if the SMR and Tempo files have the same number of total trials
% and correct trials.  If this is the case, we will assume the trial was
% recorded correctly.  ***This may change later***
htbLog = [infoStruct.tFile(1:length(infoStruct.tFile) - 3), 'log'];
[trialInfo, protocol] = TempoTrialInfo(htbLog);

SMRnumTrials = length(find([chan32.mvals] == 4));                % Total number of SMR trials.
SMRnumGoodTrials = length(find([chan32.mvals] == 12));           % Number of good trials.
tempoNumGoodTrials = length(find([trialInfo.outcomeCode] == 0 | [trialInfo.outcomeCode] == 5)); % Number of good trials in the Tempo log.

% Open the general error log so we can dump information about the
% errors for later analysis.
[x, uname] = dos('set USERNAME'); uname = uname(find(uname == '=') + 1:length(uname));
fprintf('-------------------------------------------------------------------\n');
fprintf('File: %s\n', infoStruct.smrFile);
fprintf('Protocol: %s\n', protocol);
fprintf('Date Analyzed: %s\n', date);
[x, uname] = dos('set USERNAME'); uname = uname(find(uname == '=') + 1:length(uname));
fprintf('User: %s\n', uname);

% If Tempo and Spike2 agree on the number of good trials,
% then we'll assume their were no missing marker codes.  We don't care about bad trials so if
% they're messed it is inconsequential to our data.  Otherwise, let's
% get the fixin' on.
if (tempoNumGoodTrials == SMRnumGoodTrials & length(trialInfo) == SMRnumTrials)
    fprintf('No file problems\n');
elseif (abs(length(trialInfo) - SMRnumTrials) > 5)
    fprintf('Too many total trial differences, cannot fix.\n');
    fprintf('SMR Total Trial Count: %d\nTempo Total Trial Count: %d\n', SMRnumTrials, length(trialInfo));
    fprintf('SMR Good Trial Count: %d\nTempo Good Trial Count: %d\n', SMRnumGoodTrials, tempoNumGoodTrials);
    fixable = 0;
elseif (abs(tempoNumGoodTrials - SMRnumGoodTrials) > 5)
    fprintf('Too many good trial differences, cannot fix!\n');
    fprintf('SMR Total Trial Count: %d\nTempo Total Trial Count: %d\n', SMRnumTrials, length(trialInfo));
    fprintf('SMR Good Trial Count: %d\nTempo Good Trial Count: %d\n', SMRnumGoodTrials, tempoNumGoodTrials);
    fixable = 0;
else % Fix errors
    fprintf('File errors found, attempting to fix...\n');
    smrInfo = [];
    smrIndex = 1;
    isRecovered = 1;
    
    % Finds the indices of all good trials in the event channel.  This
    % does, however, have issues in a few special cases, so there may be
    % one screwed up good trial.
    for (i = 1:length(chan32))
        % If we find a stimuls start code, then we look to see if we find
        % a success code before we get to another start code.
        if(chan32(i).mvals == 4)
            j = i + 1;
            while (j <= length(chan32) & chan32(j).mvals ~= 4)
                if (chan32(j).mvals == 12)
                    index = [index, i];
                    break;
                else
                    j = j + 1;
                end
            end
            i = j;
        end
    end
    
    % Calculate the time between the 3 and 4 codes so that when we find the
    % missing 4 trial, we can insert a 4 in the correct location.
    fprintf('Calculating average time between 3 and 4 codes for good trials\n');
    for (i = 1:length(index))
        timediff(i) = (chan32(index(i)).mark - chan32(index(i) - 1).mark) / infoStruct.chand / 25;
    end
    
    % Finds the average time between codes 3 & 4, then find the min and max
    % times.
    minMark = min(timediff); maxMark = max(timediff);
    markAve = (sum(timediff) - minMark - maxMark) / (length(timediff) - 2);
    
    % Set a flag if the min and max times vary by more than 2 milliseconds.
    if (abs(minMark - maxMark) > 2)
        fprintf('Min and Max time between codes 3 & 4 varies more than 2 milliseconds.\n');
    end
    fprintf('Average 3-4 Time: %s\n', num2str(markAve));
    
    % Separate the event channel into individual trials.
    fprintf('Separating SMR Channels...\n');
    for (i = 1:length(chan32))
        if (chan32(i).mvals == 4)
            smrInfo(smrIndex).codes = 4;
            j = i + 1;
            while (j <= length(chan32) & chan32(j).mvals ~= 4)
                smrInfo(smrIndex).codes = [smrInfo(smrIndex).codes, chan32(j).mvals];
                j = j + 1;
            end
            smrIndex = smrIndex + 1;
            i = j;
        end
    end
    
    % If the SMR trial count is 1 less than the Tempo trial count, then we
    % know we have dropped a 4 code somewhere.  The strategy to find the
    % lost 4 is to look at trials that have at least one success code and
    % check to see if they are missing a 4.
    if (length(trialInfo) ~= SMRnumTrials)
        fprintf('%d Dropped 4 codes.\n', length(trialInfo) - SMRnumTrials);
        
        isFourFound = 0;
        i = 1;
        while (i)
            % Flag any non-monotonic lines.
            %             if (length(find(diff(smrInfo(i).codes) < 0)) > 1)
            %                 fprintf(fid, 'Ambiguous trial %d: ', i); fprintf(fid, '%d ', smrInfo(i).codes); fprintf(fid, '\n');
            %             end
            
            maxEnd = max(find(smrInfo(i).codes == 8 | smrInfo(i).codes == 9 | smrInfo(i).codes == 12 | smrInfo(i).codes == 13));
            if (length(find(diff(smrInfo(i).codes(1:maxEnd)) < 0)))
                fprintf('Missing 4 code found at trial %d\n', i);
                
                if (~isempty(strfind(num2str(smrInfo(i).codes), '  15  13  ')))
                    fprintf('Found 15, 13 code combination\n');
                end
                
                % Find the 1st 3 code in the screwed up trial.  We'll put
                % the new 4 code after the 3.
                a = find(smrInfo(i).codes == 3);
                b = find([chan32.mvals] == 4);
                x.mark = markAve * 25 * infoStruct.chand + chan32(b(i) + a(1) - 1).mark;
                x.mvals = 4;
                chan32 = [chan32(1:b(i) + a(1) - 1), x, chan32(b(i) + a(1):length(chan32))];
                
                % Split the trial in smrInfo.
                y.codes = smrInfo(i).codes(1:a);
                z.codes = [4, smrInfo(i).codes(a + 1 : length(smrInfo(i).codes))];
                smrInfo = [smrInfo(1:i-1), y, z, smrInfo(i + 1:length(smrInfo))];
                
                isFourFound = 1;
            end
            i = i + 1;
            
            % This breaks us out of the loop whenever we've analyzed every
            % trial in smrInfo.
            if (i > length(smrInfo))
                i = 0;
            end
        end
        % end
        
        if (isFourFound == 0)
            fprintf('Could not find missing 4 codes\n');
        end
    else
        fprintf('No dropped 4 codes\n');
    end % End of if (length(trialInfo) ~= SMRnumTrials)
    
    % If the SMR good trial count is 1 less than the Tempo good trial
    % count, then we are missing a 12.  Finding a dropped 12 is done by
    % looking for a trial without a 12, but with a 8, 9, or 13 code.
    if (length(find([chan32.mvals] == 12)) ~= tempoNumGoodTrials)
        fprintf('%d Dropped 12 codes.\n', tempoNumGoodTrials - length(find([chan32.mvals] == 12)));
        isTwelveFound = 0;
        
        for (i = 1:length(smrInfo))
            % If we find at least one 8, 9, or 13, but no 12's in the same
            % trial, we know we've found a missing 12.
            if (length(find(smrInfo(i).codes == 8 | smrInfo(i).codes == 9 | smrInfo(i).codes == 13)) >= 1 & length(find(smrInfo(i).codes == 12)) == 0)
                % If we find a weird 15, 13 combo, we'll assume that the
                % trials sucks.  This will probably change soon.
                if (~isempty(strfind(num2str(smrInfo(i).codes), '  15  13  ')) & length(find(smrInfo(i).codes == 8 | smrInfo(i).codes == 9 | smrInfo(i).codes == 13)) == 1)
                    fprintf('Found 15, 13 combination at trial %d, with codes: ', i);
                    fprintf('%d ', smrInfo(i).codes); fprintf(' -- skipping trial...\n');
                    break;
                end
                
                fprintf('Missing 12 code found at trial %d\n', i);
                fprintf('Pre-Fix codes: '); fprintf('%d ', smrInfo(i).codes); fprintf('\n');
                
                % Because we don't really care, at this time, where the 12
                % is at, we can just turn the 8, 9, or 13 into a 12.
                b = find([chan32.mvals] == 4);
                if (length(find(smrInfo(i).codes == 13)) == 1)
                    smrInfo(i).codes(find(smrInfo(i).codes == 13)) = 12;                 
                    chan32(b(i) + find(smrInfo(i).codes == 12) - 1).mvals = 12;
                elseif (length(find(smrInfo(i).codes == 9)) == 1)
                    smrInfo(i).codes(find(smrInfo(i).codes == 9)) = 12;         
                    chan32(b(i) + find(smrInfo(i).codes == 12) - 1).mvals = 12;
                else
                    smrInfo(i).codes(find(smrInfo(i).codes == 8)) = 12;
                    chan32(b(i) + find(smrInfo(i).codes == 12) - 1).mvals = 12;
                end
                
                fprintf('Post-Fix codes: '); fprintf('%d ', smrInfo(i).codes); fprintf('\n');
                isTwelveFound = 1;
            end
        end
        
        if (isTwelveFound == 0)
            fprintf('Could not find missing 12 code\n');
            fixable = 0;
            isRecovered = 0;
        end
    else
        fprintf('No dropped 12 codes\n');
    end % End of if (length(find([CHAN32.mvals] == 12)) ~= tempoNumGoodTrials)
    
    % If we've fixed the file, we do a side by side comparison of the
    % SMR and temp data.
    if (isRecovered)
        fprintf('Performing Integrity Check...\n');
        
        % Look for all the 4/12 pairings to make sure that their total
        % matches the total of good trials found in Tempo.
        index = [];
        for (i = 1:length(chan32))
            if(chan32(i).mvals == 4)
                j = i + 1;
                while (j <= length(chan32) & chan32(j).mvals ~= 4)
                    if (chan32(j).mvals == 12)
                        index = [index, i];
                        break;
                    else
                        j = j + 1;
                    end
                end
                i = j;
            end
        end
        
        if (length(index) ~= tempoNumGoodTrials)
            fixable = 0;
            isRecovered = 0;
        end
        
        if (fixable)
            fprintf('File fixed\n');
        else
            fprintf('ERROR -- File not fixed!\n');
        end
    else
        fprintf('ERROR -- File not fixed!\n');
    end % End of if (fixed)
    
end % End of else Fix errors

% Close the file.
fprintf('-------------------------------------------------------------------\n\n\n');

return;
