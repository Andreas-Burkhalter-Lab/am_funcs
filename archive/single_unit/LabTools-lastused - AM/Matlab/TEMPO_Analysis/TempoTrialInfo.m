%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function trialInfo = TempoTrialInfo
%
%   @Author:    Christopher Broussard
%   @Date:      July, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trialInfo, protocol] = TempoTrialInfo(logFile)
trialInfo = [];
action = 'FindTrial';
numtrial = 1;

% This reads the entire contents of the Tempo log file.
logContents = textread(logFile, '%s', 'delimiter', '\n')';

% Store the protocol type.
protocol = logContents{1};

% Find every trial and records its outcome.
for (i = 1:length(logContents))
    line = strfind(logContents{i}, 'TRIAL#');
    switch action
    % Look for a Trial to analyze.
    case 'FindTrial'
        if (~isempty(line))
            action = 'FindTrialOutcome';
        end
        
    % Look for the OUTCOME line and extract the OUTCOME code.
    case 'FindTrialOutcome'
        if (isempty(line))
            line = strfind(logContents{i}, 'OUTCOME');
            if (~isempty(line))
                % Find the indexes where the outcome code exists.
                findex = strfind(logContents{i}, ' '); findex = findex(1) + 1;
                endex = strfind(logContents{i}, '('); endex = endex(1) - 1;
                
                trialInfo(numtrial).outcomeCode = str2num(logContents{i}(findex : endex));
                action = 'FindTrial';
                numtrial = numtrial + 1;
            end
        end
    end
end

return;
