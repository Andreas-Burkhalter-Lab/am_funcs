%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FindControlTrials(fileName)
%       Finds the times that are necessary to calculate the spontaneous level from a
%       Tempo log file.
%
%   @Author: Christopher Broussard
%   @Date:   April 8, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [validTrials, numtrial] = FindControlTrials(fileName)

numtrial = 0;
validTrials = [];
action = 'FindTrial';

logContents = textread(fileName, '%s', 'delimiter', '\n')';

for (i = 1:length(logContents))
    line = strfind(logContents{i}, 'TRIAL#');
    switch action    
    case 'FindTrial'      
        if (~isempty(line))
            action = 'FindCorrectTrial';
        end
    case 'FindCorrectTrial'
        if (isempty(line))
            line = strfind(logContents{i}, 'OUTCOME 0');
            if (~isempty(line))
                action = 'FindControlTrial';
                numtrial = numtrial + 1;
            end
        end
    case 'FindControlTrial'
        if (isempty(line))
            line = strfind(logContents{i}, '-9999');
            if (~isempty(line))
                action = 'FindTrial';
                validTrials = [validTrials, numtrial];
            end
        else
            action = 'FindCorrectTrial';
        end
    end
end

return;
