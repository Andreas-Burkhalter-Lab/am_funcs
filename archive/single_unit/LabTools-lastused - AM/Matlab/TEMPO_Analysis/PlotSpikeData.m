%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PlotSpikeData(handles)
%       Plots the raw spike data along with the sorted spike and window discriminator points.
%       Returns the viewable start and endpoint of the plot and handles to the threshold lines.
%
%   @Author: Christopher Broussard
%   @Date:   April 8, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [startPoint, endPoint, upperLine, lowerLine] = PlotSpikeData(handles)    
    sampleNum = 200000;     % Max number of points to be plotted.
    global CHAN1 CHAN2 CHAN32 spsData;
    
    % Determine the size of the viewing area of the plot.
    viewingArea = round(length(CHAN1) / handles.currentZoom);
    if (mod(viewingArea, 2) ~= 0)
        viewingArea = viewingArea - 1;
    end
    % Determine the domain of the midpoint min and max.
    midMin = viewingArea / 2 + 1;
    midMax = length(CHAN1) - (viewingArea / 2);
    
    % Determine the starting and ending points to plot.
    slope = midMax - midMin;
    intercept = midMin;
    x = get(handles.PlotSlider, 'Value');
    midpoint = round(slope * x + intercept);
    startPoint = midpoint - (viewingArea / 2);
    endPoint = midpoint + (viewingArea / 2);
    
    % Make sure that the endpoint and startpoint don't exceed the CHAN1 index
    % range due to rounding errors in my calculations.
    if (endPoint > length(CHAN1))
        fprintf('Shrinking endpoint from %d to %d\n', endPoint, length(CHAN1));
        endPoint = length(CHAN1);
    end
    if (startPoint < 1)
        fprintf('Increasing startpoint from %d to %d\n', startPoint, 1);
        startPoint = 1;
    end
    
    % Determine the sampling interval.
    samplingInterval = ceil(viewingArea / sampleNum);
    
    % Plot the ADC stuff
    plot([startPoint:samplingInterval:endPoint], double(CHAN1(startPoint:samplingInterval:endPoint)) / (2^15 - 1) * 5);
    axis([startPoint endPoint -handles.rangeMax handles.rangeMax]);
    hold on;

    % Plot the Event Marker stuff
    mTimes = round([CHAN32(1:length(CHAN32)).mark] / handles.chand);
    validTimes = find(mTimes >= startPoint & mTimes <= endPoint);
    if (length(validTimes) > 0)
        mHeight = [1:length(validTimes)] * 0 + handles.rangeMax - .15;
        mTimes(validTimes);
        plot(mTimes(validTimes), mHeight, 'r.');
        for(i = 1:length(validTimes))
            num = CHAN32(validTimes(i)).mvals;
            l{i} = [' ', num2str(num)];
        end
        text(mTimes(validTimes), mHeight, l);
    end
    
    % Plot the Events channel.
    var = handles.maxTime2;
    if (var > 0)
        eTimes = round([CHAN2(1:length(CHAN2)).time] / handles.chand);
        validTimes = find(eTimes >= startPoint & eTimes <= endPoint);
        eHeight = zeros(1, length(validTimes)) + handles.rangeMax - .5;
        plot(eTimes(validTimes), eHeight, '.');
    end
    
    % Plot the threshold bars.
    ubound = get(handles.UpperboundSlider, 'Value');
    lbound = get(handles.LowerboundSlider, 'Value');
    upperLine = plot([startPoint, endPoint], [ubound, ubound], 'r');
    lowerLine = plot([startPoint, endPoint], [lbound, lbound], 'm');
    
    % Plot the threshold points if necessary
    if (handles.plotThreshSpikes)
        pSpikePoints = [];
        nSpikePoints = [];
        for (i = 1:length(spsData(handles.spikeSet).spikeInfo))
            pSpikePointsIndex = find(spsData(handles.spikeSet).spikeInfo(i).pSpikeTimes > startPoint & spsData(handles.spikeSet).spikeInfo(i).pSpikeTimes < endPoint);
            nSpikePointsIndex = find(spsData(handles.spikeSet).spikeInfo(i).nSpikeTimes > startPoint & spsData(handles.spikeSet).spikeInfo(i).nSpikeTimes < endPoint);
            if (~isempty(pSpikePointsIndex))
                pSpikePoints = [pSpikePoints, spsData(handles.spikeSet).spikeInfo(i).pSpikeTimes(pSpikePointsIndex)];
            end
            if (~isempty(nSpikePointsIndex))
                nSpikePoints = [nSpikePoints, spsData(handles.spikeSet).spikeInfo(i).nSpikeTimes(nSpikePointsIndex)];
            end
        end
        
        % Plot the spikes if any exist.
        if (~isempty(pSpikePoints))
            plot(pSpikePoints, zeros(1, length(pSpikePoints)) + handles.rangeMax - .35, 'g.');
        end
        if (~isempty(nSpikePoints))
            plot(nSpikePoints, zeros(1, length(nSpikePoints)) + -handles.rangeMax + .25, 'r.');
        end
    end

    % Modify the ticks, so they represent the time in seconds rather than ADC units.
    axisTicks = get(gca, 'XTick');
    set(gca, 'XTickLabel', axisTicks / 25000.0);
    xlabel('Time in seconds'); 
    
    hold off;
    
    return;
