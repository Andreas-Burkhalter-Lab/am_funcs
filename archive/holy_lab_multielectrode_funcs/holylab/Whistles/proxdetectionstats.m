function [firstdetections,detectdurations] = proxdetectionstats (filenames)
% PROXDETECTIONSTATS: Computes statistics from the proximity detection data
%
% This function computes the time (in seconds) when the mouse first smelled the
% stimulus and the total time (in seconds) that the mouse spent smelling the
% stimulus for each *.det file.
%
% SYNTAX:
% [firstdetections,detectdurations] = proxdetectionstats (filenames)
% where
%    filenames is a cell array of *.det files
% and
%    firstdetections is a vector containing the time (in seconds) when
%      the mouse first smelled the stimulus for each file
%    detectdurations is a vector containing the total time (in seconds) 
%      that the mouse spent smelling the stimulus for each file

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>

if ~iscell(filenames)
  filenames = {filenames};
end
nfiles = length(filenames);
firstDetections = [];
detectDurations = [];
for i = 1:nfiles
    [fid,message] = fopen(filenames{i},'r');
    if (fid == -1)
        error(message);
    end
    [transitions,initVoltage,header] = ReadSensor(fid);
    
    if (header.numTransitions > 0 & initVoltage < (min([2.5 header.voltageMax])))
        secsperTransition = diff(transitions/header.scanrate);
        totaldetectTime = 0;
        if (mod(header.numTransitions,2) == 0)
            for (i=1:2:length(secsperTransition))
                totaldetectTime = totaldetectTime + secsperTransition(i);
            end
            for (i=1:2:header.numTransitions)
                if (secsperTransition(i) > 0.02)
                    firstdetectionTime = transitions(i)/header.scanrate;
                break
                else
                    firstdetectionTime = NaN;
                end
            end
            if (isnan(firstdetectionTime))
                totaldetectTime = 0;
            end
        
        else
            for (i=1:2:header.numTransitions-1)
                totaldetectTime = totaldetectTime + secsperTransition(i);
            end
        tacq = header.nscans/header.scanrate;
        lastDetection = tacq - transitions(end)/header.scanrate;
        totaldetectTime = totaldetectTime + lastDetection;
            for (i=1:2:header.numTransitions)
                if (secsperTransition(i) > 0.02)
                    firstdetectionTime = transitions(i)/header.scanrate;
                break
                else
                    firstdetectionTime = NaN;
                end
             end
            if (isnan(firstdetectionTime))
                totaldetectTime = 0;
            end
        end
    elseif (header.numTransitions == 0)
        firstdetectionTime = NaN;
        totaldetectTime = 0;
    else
        errordlg(sprintf('Proximity detection began in a high state for file %s on day %s.  You must calculate the detection time manually.',filenames{i},header.date),'Proximity Detection Error');
        firstdetectionTime = NaN;                                                                                                                 
        totaldetectTime = NaN;
    end
    firstDetections = [firstDetections;firstdetectionTime];
    detectDurations = [detectDurations;totaldetectTime];
end

firstdetections = firstDetections;
detectdurations = detectDurations;

