function DirectionRevCorrContourPlot(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;

plotWindow = figure;
set(plotWindow, 'NumberTitle', 'off');
set(plotWindow, 'Name', ['Offline Revcorr Contour Plot -- ', FILE]);
title(FILE);

g_spikeTracker = zeros(30, 30, 30, 360);
g_dirCounter = zeros(30, 30, 360);

% Set the correlation delay parameters.
corrDelayLow = 0;
corrDelayInc = 10;
corrDelayHigh = 200;

for (trial = 1:size(data.revcorr_params, 2))
    % Draw the contour subplot.
    figure(plotWindow);    
    
    fprintf('Analyzing Trial: %d\n', trial);
    
    % Grab the random directions list.
    random_directions = squeeze(data.revcorr_params(DOT_DIREC, trial, :));
    
    % Determine the unique directions from the random direction list.
    x = find(isnan(random_directions) == 0);
    uniqueDirs = munique(random_directions(1:length(x)))';
    
    % Grab the number of width and height patches.
    numXpatches = data.revcorr_params(PATCH_DIMS, 1, 1);
    numYpatches = data.revcorr_params(PATCH_DIMS, 1, 2);
    
    % Determine total number of directions per patch.
    numDirections = length(x) / (numXpatches * numYpatches);
    
    % Create the correlation delay matrix.
    corrDelay = [corrDelayLow:corrDelayInc:corrDelayHigh];
    
    % Store all the different spike histograms we collect in this data
    % structure.
    spikeHist = {};
    spikeHistIndex = 0;
    
    maxCount = 0;
    
    % Pull out the original spikes.
    origSpikes = squeeze(data.spike_data(1, :, trial));
    %origSpikes = spike_data{curr_condition,num_reps{curr_condition}}(:,1);
    
    % Figure out what the edges will be for the spike histogram.
    histEdges = find(squeeze(data.spike_data(2, :, trial)) ~= 0);
    meanEdge = round(mean(diff(histEdges)));
    histEdges = [histEdges, histEdges(length(histEdges)) + meanEdge];
    
    % Do a spike histogram for all the correlation delays.
    for (delay = corrDelay)
        spikeHistIndex = spikeHistIndex + 1;
        
        % Shift the spikes over to compensate for neuron lag.  Pad the end with
        % zeros because shifting will truncate the vector.
        shiftedSpikes = origSpikes(delay+1:length(origSpikes));
        shiftedSpikes = [shiftedSpikes, zeros(1, 5000-length(shiftedSpikes))];
        
        % Find all shiftedSpikes elements that aren't zero.
        spikeTimes = find(shiftedSpikes);
        
        % Do a histogram to bin out spike counts for each direction presented
        % in the trial.
        tmpHist = histc(spikeTimes, histEdges);
        spikeHist{spikeHistIndex} = tmpHist(1:length(tmpHist)-1);
    end
    
    % OK.  The basic strategy here is to iterate through every
    % patch and use the histogram we did above to add on directions
    % to make another histogram, which will be our analysis.
    % Example, we'll look at patch (0,0) and figure out which
    % directions were displayed in that patch during the trial.
    % Then, we'll take each direction and add it to a vector the number of
    % times it is represented in the histogram calculated above.
    for (i = 1:numXpatches)
        for (j = 1:numYpatches)           
            % Determine total number of directions per patch.
            %numDirections = totalDirections / (numXpatches * numYpatches);
            
            % Keep track of how many times each direction was displayed per
            % patch so we can normalize the results.
            for (k = 1:numDirections)
                angle = random_directions((i-1)*numDirections*numYpatches + (j-1)*numDirections + k) + 1;
                g_dirCounter(i, j, angle) = g_dirCounter(i, j, angle) + 1;
            end
            
            for (z = 1:spikeHistIndex)
                % Pull out a spike histogram from the list.
                spikeHistArray = [spikeHist{z}];
                
                for (k = 1:numDirections)
                    % Extract the direction for this segment.
                    angle = random_directions((i-1)*numDirections*numYpatches + (j-1)*numDirections + k) + 1;
                    
                    % Add the hist count for this patch, angle, delay
                    % combination.
                    g_spikeTracker(i, j, z, angle) = g_spikeTracker(i, j, z, angle) + spikeHistArray(k);
                end % End for (k = 1:numDirections)
            end % End for (z = 1:spikeHistIndex)
            
            % This keeps any of the drawing stuff from running until we're on the
            % very last trial.  This saves CPU cycles.
            if (trial == size(data.revcorr_params, 2))
                % This initializes the 2d array we'll use to make the contour
                % plot.  I add 2 to give the contour a border.
                cplot = zeros(spikeHistIndex+2, length(uniqueDirs)+2);
                
                % Loop through this patches list of correlation delayed spike
                % counts and store them in cplot to create a nice contour plot.
                for (cindex = 1:spikeHistIndex)
                    tmp = squeeze(g_spikeTracker(i, j, cindex, uniqueDirs+1))';
                    
                    % Normalize the histogram data.
                    normalizedHist = [];
                    for (div = 1:length(uniqueDirs))
                        % This makes sure we don't get a divide by zero error.
                        if (g_dirCounter(i, j, uniqueDirs(div)+1) == 0)
                            normalizedHist(div) = 0;
                        else
                            normalizedHist(div) = tmp(div) / g_dirCounter(i, j, uniqueDirs(div)+1);
                        end
                    end
                    
                    %cplot(cindex+1, :) = [0, tmp(uniqueDirs+1) / g_dirCounter(uniqueDirs+1), 0];
                    cplot(cindex+1, :) = [0, normalizedHist, 0];
                    
                    % Find the max count so we can scale all the subplots.
                    tmpMax = max(normalizedHist);
                    if (tmpMax > maxCount)
                        maxCount = tmpMax;
                    end
                end
                
                subplot(numXpatches, numYpatches, i + (j-1)*numXpatches);  
                contourf(cplot);
            end
        end % End for (j = 1:numYpatches)
    end % End for (i = 1:numXpatches)
    
    % This keeps any of the drawing stuff from running until we're on the
    % very last trial.  This saves CPU cycles.
    if (trial == size(data.revcorr_params, 2))
        % This scales all the subplots so that they are sized relative to each
        % other.
        for (i = 1:(numXpatches*numYpatches))
            subplot(numXpatches, numYpatches, i);
            %axis([-15 360 0 maxCount]);
            caxis([0, maxCount]);
            
            % Set the Y axis tick marks to represent the correlation delay.
            ticks = corrDelay(1:2:length(corrDelay));
            set(gca,'YTickLabel', ticks);
            y = [2:2:size(cplot, 1)];
            set(gca, 'YTick', y);
            set(gca, 'TickDir', 'out');
            
            % Set the X axis tick marks and labels.
            x = mean(diff(uniqueDirs));
            set(gca, 'XTick', [2:2:length(uniqueDirs)]);
            set(gca, 'XTickLabel', uniqueDirs([1:2:length(uniqueDirs)]));
            
            if (i == numXpatches*numYpatches)
                colorbar;
            elseif (i == 1)
               % Puts a label that says the filename above all the subplots.
               text(0, 30, [PATH, FILE]); 
            end
        end
    end
end % End for (trial = 1:size(data.revcorr_params, 2))

return;