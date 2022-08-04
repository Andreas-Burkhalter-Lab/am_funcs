function [flagged_trials trialdata_out nflaggedTrials] = manual_trial_noise_check(trialdata,snipdata,opt_trialCheck,merec_obj,window_scans)
%MANUAL_TRIAL_NOISE_CHECK Manually remove noisy high spiking trials.
%  Traces from high-spiking trials and channels will be presented as traces; user can accept
%  or reject each trial. Ouput 'noisy trials' gives the indices of trials
%  from trialdata marked by user for rejection.
%%%     options struct:
%%% opt_trialCheck.zScoreThresh.... above-average-spiking trials with this z-score or greater will be flagged
%%% opt_trialCheck.nChansThresh... trials with at least this many flagged channels will be presented for checking
%%% opt_trialCheck.maxRows... max number of subplot rows to use for plotting traces - fill all rows before columns
%%% opt_trialCheck.maxColumns... max number of subplot rows to use for plotting traces
%%% opt_trialCheck.chans... channels to check; defaults to all channels if empty
%%% opt_trialCheck.plotTrialWaveforms... plot the full waveform of the trial alongside the snips
%
%%% merec_obj and window_scans only used if opt_trialCheck.plotTrialWaveforms == 1; merec
%%%     file will be taken from snipdata filehandle if not specified
% last updated 10/23/15 on vivid

if (~exist('merec_obj','var') || isempty(merec_obj)) && opt_trialcheck.plotTrialWaveforms 
    merec_obj = merecmm(fullfile(snipdata.fh.abspathstr,snipdata.fh.filename));
end

%% Mark and load potentially noisy trials on channels of interest.
if isempty(opt_trialCheck.chans) % default to all channels if unspecified
    opt_trialCheck.chans = snipdata.channels; 
end
chanIndices = NaN(1,length(opt_trialCheck.chans));
for inputChanInd = 1:length(chanIndices) % get indices of channels specified in opt_trialCheck.chans
    chanIndices(inputChanInd) = find(opt_trialCheck.chans(inputChanInd) == snipdata.channels);
end
    
ntrials = size(trialdata,1);
flagged_trials = NaN(ntrials,1);
trialdata.do_check = false(ntrials,1);
trialdata.flagged = NaN(ntrials,1); % 1 = noisy, 0 = not noisy (checked), NaN = not noisy (not checked)
trialdata.nHighSpikeChs = NaN(ntrials,1);
trialdata.highSpikeChInds = cell(ntrials,1); 
trialdata.zscore_spikes = NaN(size(trialdata.spikes)); % keep unused channels as NaNs so that channel indices remain valid
trialdata.zscore_spikes(:,chanIndices) = zscore(trialdata.spikes(:,chanIndices));
for indTrial = 1:ntrials     % list all channels exceeding the zscore threshold
    trialdata.highSpikeChInds{indTrial} = find(trialdata.zscore_spikes(indTrial,:) >= opt_trialCheck.zScoreThresh); % channel INDICES, not names
end
trialdata.nHighSpikeChs = cell2mat(cellfun(@length,trialdata.highSpikeChInds,'UniformOutput',false)); % number of high spiking channels
trialdata.do_check = trialdata.nHighSpikeChs > opt_trialCheck.nChansThresh; % logicals
nTrialsToCheck = length(find(trialdata.do_check));
trialsToCheck = find(trialdata.do_check); % trial numbers to check
chanIndsToCheck = unique(cell2mat(trialdata.highSpikeChInds'));  % all channel INDICES, not names, to check from all trials
nChansToCheck = length(chanIndsToCheck); % not all of these may be checked if there are more than maxRows*maxColumns

% Load (but don't display) snips from all channels for trials of interest
% to keep channel indices valid. 
snipsToCheck = snips_from_trialdata(trialdata,snipdata,trialsToCheck,opt_trialCheck.chans); % dataset... maybe don't need to get all of these snips
fprintf('Found %g/%g trials to check (z-score thresh = %g, nChans thresh = %g).\n',...
    nTrialsToCheck, ntrials, opt_trialCheck.zScoreThresh, opt_trialCheck.nChansThresh)

%% Manual trial checking
% Keep the same subplot layout from trial to trial so that the same channel
% is always plotted in the same location at the same size.
if ~isfield(snipdata,'unitnames') % add unitnames field to snipdata if not present
    snipdata.unitnames = cellfun(@(x)['ch',num2str(x)],num2cell(snipdata.channels),'UniformOutput',false);
elseif isfield(snipdata,'unitnames') && length(snipdata.channels) ~=length(snipdata.unitnames)       %%%%%%% why is this check necessary? 
    error('manual_trial_noise_check.m assumes that length(snipdata.channels)==length(snipdata.unitnames); here they are unequal.')
end
nChansToPlotPerTrial = min([nChansToCheck, opt_trialCheck.maxRows*opt_trialCheck.maxColumns]);  
chanIndsToPlot = chanIndsToCheck(1:nChansToPlotPerTrial); % INDICES, not names; plot and check no more than nChansToPlot channels
plotRows = min([nChansToPlotPerTrial, opt_trialCheck.maxRows]); % number of subplot rows; make all rows before making a second column
plotColumns = min([ceil(nChansToPlotPerTrial/plotRows), opt_trialCheck.maxColumns]); % number of subplot columns
%%% maintain fixed positions for channels within the plot there are at
%%% least as many plot spaces as channels to plot
if nChansToCheck <= plotRows*plotColumns
    fixedChanPlotPositions = 1;
else
    fixedChanPlotPositions = 0;
end

if opt_trialCheck.plotTrialWaveforms
    plotColumns = 2*plotColumns; % add extra columns for the full trial waveforms
end
origFigWindowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','docked'); % dock the figure so it's not covered over if commandwindow is called
hfig = figure;

indCheckTrial = 0; % counter for indices within snipsToCheck, not trial numbers from trialdata
while indCheckTrial <= nTrialsToCheck-1
   indCheckTrial = indCheckTrial + 1;
   thisTrial = trialsToCheck(indCheckTrial); % trial number from trialdata
   clf(hfig);
   nChansToPlotThisTrial = min([nChansToPlotPerTrial length(trialdata.highSpikeChInds{thisTrial})]); 
   chansToPlotThisTrial = snipdata.channels(trialdata.highSpikeChInds{thisTrial}(1:nChansToPlotThisTrial));
   for chanName = chansToPlotThisTrial % % 
       indChan = find(chanName == snipdata.channels); 
       tempSnips = snipsToCheck{indCheckTrial,[snipdata.unitnames{indChan} '_snips']}; % snips to plotclc
       
       nsnips = size(tempSnips,2);
       temptitle = [snipdata.unitnames{indChan} ' ... ' num2str(nsnips) ' snips.'];
       if fixedChanPlotPositions
           plotIndex = find(indChan==chanIndsToPlot); % index within chanIndsToPlot... maintain location for each chan
       else
           plotIndex = find(chanName == chansToPlotThisTrial); % location depends on index within chs to plot this trial
       end
       if opt_trialCheck.plotTrialWaveforms
            trialWave = merec_obj(chanName,[trialdata.onset(thisTrial):trialdata.onset(thisTrial)+window_scans]); 
            thisRow = rem(plotIndex-1,plotRows)+1;
            thisColumn = 2*(floor([plotIndex-1]/plotRows) + 1) - 1;
            plotLocSnips = sub2ind([plotColumns plotRows], thisColumn, thisRow); % alternate snips and waves columns
            plotLocTrialwaves = plotLocSnips + 1; % put waves to the right of corresponding snips
            subplot(plotRows,plotColumns,plotLocTrialwaves); 
            plot(trialWave)
            title(temptitle);
       else
            plotLocSnips = plotIndex; % plotting snips only, not full trial waveforms
       end
       subplot(plotRows,plotColumns,plotLocSnips);
       plot(tempSnips);
       title(temptitle)
   end
   flagAnswer = input(sprintf(['Trial %g (%g/%g): enter ''d'' to discard , ',...
       '[blank] to keep, or # to go back # trials.\n '],...
       thisTrial, indCheckTrial, nTrialsToCheck),'s');
   numAnswer = str2num(flagAnswer);
   if strcmp(flagAnswer,'d')
       trialdata.flagged(thisTrial) = 1; % mark trial for discarding
   elseif isempty(flagAnswer)
       trialdata.flagged(thisTrial) = 0; % keep this trial, use 0 to mark that it was checked
   elseif ~isempty(numAnswer) && numAnswer>0 && rem(numAnswer,1)==0  % if a positive integer was entered
        indCheckTrial = max([0 indCheckTrial-1-numAnswer]); % go back flagAnswer trials or to first trial
   else
       fprintf('Input not recognized.\n')
       indCheckTrial = indCheckTrial - 1; % try the same trial again to get proper input
   end
end
       
flagged_trials = find(trialdata.flagged==1); % trial numbers marked as noisy
flagged_trials = flagged_trials(~isnan(flagged_trials)); % discard unused indices
fprintf('%g/%g trials discarded, % checked.\n',length(flagged_trials),ntrials,nTrialsToCheck);

if nargout >= 2
    trialdata_out = trialdata;
end
if nargout >= 3
    nflaggedTrials = length(flagged_trials);
end
   
close(hfig);
set(0,'DefaultFigureWindowStyle',origFigWindowStyle);