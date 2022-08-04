%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function fixable = FixSMRData(handles)
%       Finds and fixes problems with SMR data that disagrees with Tempo
%       data.  Returns a 1 if the SMR file is fixable, 0 if it's not.
%
%   @Author:    Christopher Broussard
%   @Date:      July, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fixable, isRecovered] = FixSMRData(handles)
global CHAN32;
fixable = 1;
isRecovered = 0;
index = [];

% Check to see if the SMR and Tempo files have the same number of total trials
% and correct trials.  If this is the case, we will assume the trial was
% recorded correctly.  ***This may change later***
htbLog = [handles.htdFile(1:length(handles.htdFile) - 3), 'log'];
[trialInfo, protocol] = TempoTrialInfo(htbLog);

SMRnumTrials = length(find([CHAN32.mvals] == 4));                % Total number of SMR trials.
SMRnumGoodTrials = length(find([CHAN32.mvals] == 12));           % Number of good trials.
tempoNumGoodTrials = length(find([trialInfo.outcomeCode] == 0 | [trialInfo.outcomeCode] == 5)); % Number of good trials in the Tempo log.

% Open the general error log so we can dump information about the
% errors for later analysis.
[x, uname] = dos('set USERNAME'); uname = uname(find(uname == '=') + 1:length(uname));
fid = fopen(['Z:/Data/SpikeSort Errors/spikesort_errors -- ', uname(1:length(uname) - 1), '.log'], 'a');
if (fid == -1)
    disp('Could not open error log');
else
    fprintf(fid, '-------------------------------------------------------------------\n');
    fprintf(fid, 'File: %s\n', ['Z:/Data/CED/', handles.monkeyName, '/', handles.dataFileName]);
    fprintf(fid, 'Protocol: %s\n', protocol);
    fprintf(fid, 'Date Analyzed: %s\n', date);
    [x, uname] = dos('set USERNAME'); uname = uname(find(uname == '=') + 1:length(uname));
    fprintf(fid, 'User: %s\n', uname);
end

% If Tempo and Spike2 agree on the number of good trials,
% then we'll assume their were no missing marker codes.  We don't care about bad trials so if
% they're messed it is inconsequential to our data.  Otherwise, let's
% get the fixin' on.
if (tempoNumGoodTrials == SMRnumGoodTrials & length(trialInfo) == SMRnumTrials)
    fprintf(fid, 'No file problems\n');
elseif (abs(length(trialInfo) - SMRnumTrials) > 5)
    fprintf(fid, 'Too many total trial differences, cannot fix.\n');
    fprintf(fid, 'SMR Total Trial Count: %d\nTempo Total Trial Count: %d\n', SMRnumTrials, length(trialInfo));
    fprintf(fid, 'SMR Good Trial Count: %d\nTempo Good Trial Count: %d\n', SMRnumGoodTrials, tempoNumGoodTrials);
    fixable = 0;
elseif (abs(tempoNumGoodTrials - SMRnumGoodTrials) > 5)
    fprintf(fid, 'Too many good trial differences, cannot fix!\n');
    fprintf(fid, 'SMR Total Trial Count: %d\nTempo Total Trial Count: %d\n', SMRnumTrials, length(trialInfo));
    fprintf(fid, 'SMR Good Trial Count: %d\nTempo Good Trial Count: %d\n', SMRnumGoodTrials, tempoNumGoodTrials);
    fixable = 0;
else % Fix errors
    fprintf(fid, 'File errors found, attempting to fix...\n');
    smrInfo = [];
    smrIndex = 1;
    isRecovered = 1;
    
    % Finds the indices of all good trials in the event channel.  This
    % does, however, have issues in a few special cases, so there may be
    % one screwed up good trial.
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
    
    % Calculate the time between the 3 and 4 codes so that when we find the
    % missing 4 trial, we can insert a 4 in the correct location.
    fprintf(fid, 'Calculating average time between 3 and 4 codes for good trials\n');
    for (i = 1:length(index))
        timediff(i) = (CHAN32(index(i)).mark - CHAN32(index(i) - 1).mark) / handles.chand / 25;
    end
    
    % Finds the average time between codes 3 & 4, then find the min and max
    % times.
    minMark = min(timediff); maxMark = max(timediff);
    markAve = (sum(timediff) - minMark - maxMark) / (length(timediff) - 2);
    
    % Set a flag if the min and max times vary by more than 2 milliseconds.
    if (abs(minMark - maxMark) > 2)
        fprintf(fid, 'Min and Max time between codes 3 & 4 varies more than 2 milliseconds.\n');
    end
    fprintf(fid, 'Average 3-4 Time: %s\n', num2str(markAve));
    
    % Separate the event channel into individual trials.
    fprintf(fid, 'Separating SMR Channels...\n');
    for (i = 1:length(CHAN32))
        if (CHAN32(i).mvals == 4)
            smrInfo(smrIndex).codes = 4;
            j = i + 1;
            while (j <= length(CHAN32) & CHAN32(j).mvals ~= 4)
                smrInfo(smrIndex).codes = [smrInfo(smrIndex).codes, CHAN32(j).mvals];
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
        fprintf(fid, '%d Dropped 4 codes.\n', length(trialInfo) - SMRnumTrials);
        
        isFourFound = 0;
        i = 1;
        while (i)
            % Flag any non-monotonic lines.
            %             if (length(find(diff(smrInfo(i).codes) < 0)) > 1)
            %                 fprintf(fid, 'Ambiguous trial %d: ', i); fprintf(fid, '%d ', smrInfo(i).codes); fprintf(fid, '\n');
            %             end
            
            maxEnd = max(find(smrInfo(i).codes == 8 | smrInfo(i).codes == 9 | smrInfo(i).codes == 12 | smrInfo(i).codes == 13));
            if (length(find(diff(smrInfo(i).codes(1:maxEnd)) < 0)))
                fprintf(fid, 'Missing 4 code found at trial %d\n', i);
                
                if (~isempty(strfind(num2str(smrInfo(i).codes), '  15  13  ')))
                    fprintf(fid, 'Found 15, 13 code combination\n');
                end
                
                % Find the 1st 3 code in the screwed up trial.  We'll put
                % the new 4 code after the 3.
                a = find(smrInfo(i).codes == 3);
                b = find([CHAN32.mvals] == 4);
                x.mark = markAve * 25 * handles.chand + CHAN32(b(i) + a(1) - 1).mark;
                x.mvals = 4;
                CHAN32 = [CHAN32(1:b(i) + a(1) - 1), x, CHAN32(b(i) + a(1):length(CHAN32))];
                
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
            fprintf(fid, 'Could not find missing 4 codes\n');
        end
    else
        fprintf(fid, 'No dropped 4 codes\n');
    end % End of if (length(trialInfo) ~= SMRnumTrials)
    
    % If the SMR good trial count is 1 less than the Tempo good trial
    % count, then we are missing a 12.  Finding a dropped 12 is done by
    % looking for a trial without a 12, but with a 8, 9, or 13 code.
    if (length(find([CHAN32.mvals] == 12)) ~= tempoNumGoodTrials)
        fprintf(fid, '%d Dropped 12 codes.\n', tempoNumGoodTrials - length(find([CHAN32.mvals] == 12)));
        isTwelveFound = 0;
        
        for (i = 1:length(smrInfo))
            % If we find at least one 8, 9, or 13, but no 12's in the same
            % trial, we know we've found a missing 12.
            if (length(find(smrInfo(i).codes == 8 | smrInfo(i).codes == 9 | smrInfo(i).codes == 13)) >= 1 & length(find(smrInfo(i).codes == 12)) == 0)
                % If we find a weird 15, 13 combo, we'll assume that the
                % trials sucks.  This will probably change soon.
                if (~isempty(strfind(num2str(smrInfo(i).codes), '  15  13  ')) & length(find(smrInfo(i).codes == 8 | smrInfo(i).codes == 9 | smrInfo(i).codes == 13)) == 1)
                    fprintf(fid, 'Found 15, 13 combination at trial %d, with codes: ', i);
                    fprintf(fid, '%d ', smrInfo(i).codes); fprintf(fid, ' -- skipping trial...\n');
                    break;
                end
                
                fprintf(fid, 'Missing 12 code found at trial %d\n', i);
                fprintf(fid, 'Pre-Fix codes: '); fprintf(fid, '%d ', smrInfo(i).codes); fprintf(fid, '\n');
                
                % Because we don't really care, at this time, where the 12
                % is at, we can just turn the 8, 9, or 13 into a 12.
                b = find([CHAN32.mvals] == 4);
                if (length(find(smrInfo(i).codes == 13)) == 1)
                    smrInfo(i).codes(find(smrInfo(i).codes == 13)) = 12;                 
                    CHAN32(b(i) + find(smrInfo(i).codes == 12) - 1).mvals = 12;
                elseif (length(find(smrInfo(i).codes == 9)) == 1)
                    smrInfo(i).codes(find(smrInfo(i).codes == 9)) = 12;         
                    CHAN32(b(i) + find(smrInfo(i).codes == 12) - 1).mvals = 12;
                else
                    smrInfo(i).codes(find(smrInfo(i).codes == 8)) = 12;
                    CHAN32(b(i) + find(smrInfo(i).codes == 12) - 1).mvals = 12;
                end
                
                fprintf(fid, 'Post-Fix codes: '); fprintf(fid, '%d ', smrInfo(i).codes); fprintf(fid, '\n');
                isTwelveFound = 1;
            end
        end
        
        if (isTwelveFound == 0)
            fprintf(fid, 'Could not find missing 12 code\n');
            fixable = 0;
            isRecovered = 0;
        end
    else
        fprintf(fid, 'No dropped 12 codes\n');
    end % End of if (length(find([CHAN32.mvals] == 12)) ~= tempoNumGoodTrials)
    
    % If we've fixed the file, we do a side by side comparison of the
    % SMR and temp data.
    if (isRecovered)
        fprintf(fid, 'Performing Integrity Check...\n');
        
        % Look for all the 4/12 pairings to make sure that their total
        % matches the total of good trials found in Tempo.
        index = [];
        for (i = 1:length(CHAN32))
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
        
        if (length(index) ~= tempoNumGoodTrials)
            fixable = 0;
            isRecovered = 0;
        end
        
        if (fixable)
            fprintf(fid, 'File fixed\n');
        else
            fprintf(fid, 'ERROR -- File not fixed!\n');
        end
    else
        fprintf(fid, 'ERROR -- File not fixed!\n');
    end % End of if (fixed)
    
end % End of else Fix errors

% Close the file.
fprintf(fid, '-------------------------------------------------------------------\n\n\n');
fclose(fid);

return;
