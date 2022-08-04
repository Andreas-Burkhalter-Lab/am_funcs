%%% get spike timing info and determine whether a cell significantly
%%% responds to the stimulus relative to baseline
%%%% called by getTempoTrialData.m after getSpikeTimes.m
% last saved 8/29/16

function analyze_spike_timing(loomingSpeedDependentSpikeWindow)
% % % % % % %%% data = data struct from LoadTEMPOData.m
%%% loomingSpeedDependentSpikeWindow: if on, spike counting window for looming stim will depend on looming speed
%
% spikeTimingFile is the psth file for all runs created by
% getTempoTrialData
postOnsetWindowStart = 50; % beginning of spike counting window after stim onset in ms   
postOnsetWindowEnd = 500; % end of spike counting window after stim onset in ms    
event_bin_width = 0.001; % duration in seconds of a time bin for data in spikesTable
signCutoff = 0.05; % p-value cutoff for determining significant difference between baseline and response
spikeTimingFile = 'C:\Users\AM\Documents\Lab\recordings\spike_timing_summary.mat';

load(spikeTimingFile,'spikes_list','column_headings');
responseWindowDuration = postOnsetWindowEnd-postOnsetWindowStart; % ms

% check for duplicates in spikes_list
fnamelist = spikes_list(:,1);
[~,unqInds] = unique(flipud(fnamelist),'stable'); % get unique names starting from the bottom to keep newest entries, discard older
unqInds = length(fnamelist) - unqInds + 1;        % get unique names starting from the bottom to keep newest entries, discard older
if length(unique(spikes_list(:,1))) < size(spikes_list,1) % if there are duplicates
    spikes_list = spikes_list(unqInds); % keep only the first of each duplicate
    save(spikeTimingFile,'spikes_list','-append'); % save over spikes_list in the file containing duplicates
end

spikesTable = cell2table(spikes_list,'VariableNames',column_headings);
nruns = height(spikesTable);
spikesTable.baseHz = NaN(nruns,1);
spikesTable.respHz = NaN(nruns,1);
spikesTable.baseSD = NaN(nruns,1);
spikesTable.respSD = NaN(nruns,1);
spikesTable.resp_p = NaN(nruns,1);
spikesTable.sgn_resp = NaN(nruns,1); % 0=nonsig response, 1 = sig response


for runInd = 1:nruns
    onset = spikesTable.stim_onset_time(runInd);
%     offset = spikesTable.stim_offset_time(runInd);
    if ~isempty(spikesTable.spikes_null) && isnumeric(spikesTable.spikes_null{1,1}) % fix formatting
        spnull = rowfun(@(x){x},spikesTable(:,'spikes_null'));
        spikesTable.spikes_null = spnull{:,:}; 
    end
    if ~isempty(spikesTable.spikes_opt) && isnumeric(spikesTable.spikes_opt{1,1}) % fix formatting
        spopt = rowfun(@(x){x},spikesTable(:,'spikes_opt'));
        spikesTable.spikes_opt = spopt{:,:}; 
    end
    
    if ~isempty(spikesTable.spikes_null)
        theseNullSpikes = spikesTable.spikes_null{runInd};
    elseif isempty(spikesTable.spikes_null)
        theseNullSpikes = {};
    end
    if ~isempty(spikesTable.spikes_opt)
        theseOptSpikes = spikesTable.spikes_opt{runInd};
    elseif isempty(spikesTable.spikes_opt)
        theseOptSpikes = {};
    end

    nNullTrials = length(theseNullSpikes);
    nOptTrials = length(theseOptSpikes);
    thisBaseline = NaN(size(theseNullSpikes));
    thisResponse = NaN(size(theseOptSpikes));   
    if spikesTable.has_null_trials(runInd) % if run has null trials, get baseline from null trials
        for trialInd = 1:nNullTrials
            null_spikes = length(find([theseNullSpikes{trialInd} >= onset+postOnsetWindowStart] & [theseNullSpikes{trialInd} <= onset+postOnsetWindowEnd]));
            thisBaseline(trialInd) = null_spikes/responseWindowDuration ./ event_bin_width; % onset is in ms, multiply by 1000 to get s (hz)
        end
    end
    for trialInd = 1:nOptTrials % get response spikes and if there aren't null trials, get baseline spikes 
        if ~spikesTable.has_null_trials(runInd) % if there aren't null trials
            base_spikes = length(find(theseOptSpikes{trialInd} <= onset));
            thisBaseline(trialInd) = base_spikes/onset ./ event_bin_width; % onset is in ms, multiply by 1000 to get seconds (hz)
        end
        % if it's a looming run and loomingSpeedDependentSpikeWindow == 1,
        % do not set new spike response windows, use the ones already set
        % by ComputeSpikeRates_looming.m in
        if loomingSpeedDependentSpikeWindow && any(strcmp(spikesTable.protocol_name(runInd),{'LOOMING_SPEED';'LOOMING_GAMP';'LOOMING_GASSIAN_AMP'})) % if a looming trial and loomingSpeedDependentSpikeWindow == 1
            thisResponse(trialInd) = spikesTable.spikes_loom_rsp{runInd}(trialInd);
        else
            resp_spikes = length(find([theseOptSpikes{trialInd} >= onset+postOnsetWindowStart] & [theseOptSpikes{trialInd} <= onset+postOnsetWindowEnd]));
            thisResponse(trialInd) = resp_spikes/responseWindowDuration * 1000; % get response spike rate in hz; onset is in ms, multiply by 1000 to get seconds (hz)
        end
    end
    spikesTable.baseHz(runInd) = mean(thisBaseline);
    spikesTable.respHz(runInd) = mean(thisResponse);
    spikesTable.baseSD(runInd) = std(thisBaseline);
    spikesTable.respSD(runInd) = std(thisResponse);
    if spikesTable.has_null_trials(runInd) % unpaired ttest if rates are from different trials
        [~,spikesTable.resp_p(runInd)] = ttest2(thisBaseline, thisResponse);        
    elseif ~spikesTable.has_null_trials(runInd) % paired ttest if rates are from the same trial
        [~,spikesTable.resp_p(runInd)] = ttest(thisBaseline, thisResponse);
    end
    spikesTable.sgn_resp(runInd) = spikesTable.resp_p(runInd) < signCutoff;
end

save(spikeTimingFile,'spikesTable','-append'); % save analysis