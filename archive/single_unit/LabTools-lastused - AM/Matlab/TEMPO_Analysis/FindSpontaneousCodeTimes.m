%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FindSpontaneousCodeTimes(handles)
%       Finds the times that are necessary to calculate the spontaneous level.
%
%   @Author: Christopher Broussard
%   @Date:   April 8, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function codeTimes = FindSpontaneousCodeTimes(handles)
global CHAN32;
index = [];
controlTrials = [];
tCount = 1;

% Find all the stimulus start codes that resulted in a successful trial.
for (i = 1:length(CHAN32))
    % If we find a stimulus start code, then we look to see if we find
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

% Check to see if the Tempo log indicates there are trials which will act
% as a control.  If the file doesn't exist, then calculate the spontaneous
% code times and anything between a 3 & 4 code.  Otherwise, the codetimes are
% based on the time between 4 & 5 codes only in control trials specified in
% the log file.
htbLog = [handles.htdFile(1:length(handles.htdFile) - 3), 'log'];
if (exist(htbLog, 'file') > 0)
    [controlTrials, junk] = FindControlTrials(htbLog);
end

if (isempty(controlTrials))
    % Find all time periods between a fix code and a stim start code.
    for (i = index - 1)
        % If we find a fixate code, record the start and stop times.
        if(CHAN32(i).mvals == 3)
            j = i + 1;
            while (j <= length(CHAN32) & CHAN32(j).mvals ~= 3)
                if (CHAN32(j).mvals == 4)
                    codeTimes(tCount).sTime = round(CHAN32(i).mark / handles.chand);
                    codeTimes(tCount).eTime = round(CHAN32(j).mark / handles.chand);
                    tCount = tCount + 1;
                    break;
                else
                    j = j + 1;
                end
            end
            i = j;
        end
    end
else
    % Find all the code times for the control trials in the tempo log.
    for (i = index(controlTrials))
        % Record the start time of a control trial.
        codeTimes(tCount).sTime = round(CHAN32(i).mark / handles.chand);
        
        % Find the following 5 code, and record its time.
        for (j = i + 1:length(CHAN32))
            if (CHAN32(j).mvals == 5)
                codeTimes(tCount).eTime = round(CHAN32(j).mark / handles.chand);
                tCount = tCount + 1;
                break;
            end
        end
    end
end

return;
