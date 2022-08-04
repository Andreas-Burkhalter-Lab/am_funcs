%RF_GETSPIKES: Get spike snippets from the electrode data for
%analyze_rf_responses.
%%% Last updated 1/30/16 on vivid

opt_autosnip2.files = {[opt_autosnip2.files '.merec']}; % needs to be cell for autosnip2 to handle it correctly

do_autosnip2_autosort = 1;
if ~rerun_autosnip2_autosort
    if exist([savedir filesep 'overview.mat'],'file')
        load([savedir filesep 'overview.mat']);
        if strcmp(getfname(sorthead.fh.filename), getfname(opt_autosnip2.files{:}))
            disp('Found autosort overview file; not rerunning autosnip2 and autosort.')
            do_autosnip2_autosort = 0;
        end
    end
end
if do_autosnip2_autosort
    autosnip2_AM(opt_autosnip2, opt_snippetfile); % opt_snippetfile specifies filtering parameters
    autosort_AM(snipfile2sortheader(strcat(opt_autosnip2.files{:}(1:end-6),'.ssnp')), savedir);
    load([savedir filesep 'overview.mat']); % load newly created overview.mat
end

snipdata = sorthead_from_raw(sorthead, rec_chans, []);

%% Count stimulus-evoked spike snippets after each stimulus event
% Find mean number of spike snippets in the response window
% following each stimulus type. All times should be measured in scans. 
%%% counts_mean is a vector of length = # unique event labels listing the 
%%%     number of spikes following  the corresponding event within time 
%%%     'window', averaged over all iterations
%%% counts represents spiking responses as a grid for each channel for each
%%%     iteration
%%% trialdata_rf represents each trial as a row in a dataset
%%%     .snip_ind lists the indices (out of all snips on this channel from .ssnp file) 
%%%         of snips occuring in this trial window 
events_per_iter = nevents/stimpars_rf.nAngles;
response_time_total = (size(trialdata_rf,1)*resp_window); % total time (s) counted as response
baseline_time_total = (merec_obj.nscans / merec_obj.scanrate) - response_time_total; % total time (s) counted as baseline
nangrid = {NaN(stimpars_rf.rows, stimpars_rf.columns)}; 
chanNames = [cellfun(@(x)['ch_',num2str(x)],num2cell(rec_chans),'UniformOutput',false) 'avg'];
iterNames = [cellfun(@(x)['iter_',num2str(x)],num2cell(1:stimpars_rf.nAngles),'UniformOutput',false) 'avg'];
counts = dataset({repmat(nangrid,length(rec_chans)+1,stimpars_rf.nAngles+1), iterNames{:}}, ...
    {NaN(length(chanNames),4),'response_hz','baseline_hz','response_total','baseline_total'}, ...
    'ObsNames', chanNames);
trialdata_rf.iter = NaN(nevents,1); % iteration of the full grid presentation
trialdata_rf.row = NaN(nevents,1); % stimulus grid row
trialdata_rf.column = NaN(nevents,1); % stimulus grid column
trialdata_rf.spikes = NaN(nevents,nchans); 
trialdata_rf.snip_ind = cell(size(trialdata_rf,1),length(snipdata.channels));

% Count spikes on each channel on each trial.
% could shift window back by ~0.5ms for first half of spike
for indChan = 1:length(rec_chans) % channel INDICES, not channel names
    tsnips = snipdata.sniptimes{indChan};
    thisChanName = chanNames(indChan);
    countsmat = NaN(stimpars_rf.rows, stimpars_rf.columns, stimpars_rf.nAngles); % initialize and clear
    for indTrial = 1:nevents  %% for each event, get number of snips following that event on this channel
        iteration = floor((indTrial-1)/(stimpars_rf.rows*stimpars_rf.columns)) + 1; 
        row = floor(event_bins(3,indTrial)/max_vsteps);
        column = mod(event_bins(3,indTrial),max_vsteps);
        trialStartScan = trialdata_rf.onset(indTrial); % scan when the trial-start pulse was detected
        trialFinishedScan = trialStartScan + window_scans; 
        trialdata_rf.iter(indTrial) = iteration;
        trialdata_rf.row(indTrial) = row; % stimulus grid row
        trialdata_rf.column(indTrial) = column; 
        trialdata_rf.snip_ind{indTrial,indChan} = find(tsnips > trialStartScan & tsnips <= trialFinishedScan);
        nSpikesThisTrial = length(trialdata_rf.snip_ind{indTrial,indChan}); % count spikes 
        trialdata_rf.spikes(indTrial,indChan) = nSpikesThisTrial; 
        if opt_trialCheck.keepRaw % save spike/trial data before removing noisy trials
            counts{chanNames(indChan),iterNames(iteration)}(row,column) = nSpikesThisTrial; % copy data to counts
        end  
    end
    if opt_trialCheck.keepRaw % only compute counts at this point if the raw data is going to be saved
        for iter = 1:stimpars_rf.nAngles
            countsmat(:,:,iter) = counts{chanNames(indChan),iterNames(iter)};  
        end
        counts.avg{thisChanName} = mean(countsmat,3); % average for this channel over all iterations
        counts.response_total(thisChanName) = sum(sum(sum(countsmat))); % total snips in response windows for this chan 
        counts.response_hz(thisChanName) = counts.response_total(thisChanName) / response_time_total; % avg hz during resp window
        counts.baseline_total(thisChanName) = size(snipdata.snips{indChan},2) - counts.response_total(thisChanName); % total snips during baseline
        counts.baseline_hz(thisChanName) = counts.baseline_total(thisChanName) / baseline_time_total; % avg hz during baseline
    end
end
trialdata_rf.rowcolumn = []; % eliminate unused variable

% Manual trial checking
if manual_trial_check % add variable for spikes automatically detected, with manually removed trials turned to NaNs
    [flagged_trials trialnoise_rf nFlaggedTrials] = manual_trial_noise_check(trialdata_rf,snipdata,opt_trialCheck,merec_obj,window_scans);
    trialdata_rf.discard = trialnoise_rf.flagged; % if =1, this trial will not be counted
    trialdata_rf = trialdata_rf(:,{'discard','iter','row','column','onset','dur','spikes','snip_ind'}); % reorder variables
    if opt_trialCheck.keepRaw % save spike/trial data before removing noisy trials
        counts_raw = counts; 
        trialdata_rf_raw = trialdata_rf;
    end
    trialdata_rf.spikes(trialdata_rf.discard==1,:) = NaN; % set spike data to NaN for all trials flagged as noisy
else
    clear opt_trialcheck
end

% Create 'counts' grid representation of spike data after (optionally) noisy
% trials have been blanked. 
for indChan = 1:length(snipdata.sniptimes) % channel INDICES, not channel names
    for indTrial = 1:nevents  %% for each event, get number of snips following that event on this channel    
        iteration = trialdata_rf.iter(indTrial);
        row = trialdata_rf.row(indTrial); % stimulus grid row
        column = trialdata_rf.column(indTrial);
        thisChanName = chanNames(indChan);
        nSpikesThisTrial = trialdata_rf.spikes(indTrial,indChan);
        counts{chanNames(indChan),iterNames(iteration)}(row,column) = nSpikesThisTrial; % copy data to counts
    end
    for iter = 1:stimpars_rf.nAngles
        countsmat(:,:,iter) = counts{chanNames(indChan),iterNames(iter)};  
    end
    counts.avg{thisChanName} = nanmean(countsmat,3); % average for this channel over all iterations
    counts.response_total(thisChanName) = nansum(sum(sum(countsmat))); % total snips in response windows for this chan 
    counts.response_hz(thisChanName) = counts.response_total(thisChanName) / response_time_total; % avg hz during resp window
    counts.baseline_total(thisChanName) = size(snipdata.snips{indChan},2) - counts.response_total(thisChanName); % total snips during baseline
    counts.baseline_hz(thisChanName) = counts.baseline_total(thisChanName) / baseline_time_total; % avg hz during baseline
end

% Get averages for a given trial across all channels and grand-average.
countsmat = NaN(stimpars_rf.rows, stimpars_rf.columns, nchans, stimpars_rf.nAngles); % initialize and clear
for iter = 1:stimpars_rf.nAngles
    for indChan = 1:length(snipdata.sniptimes)
        thisChanName = chanNames(indChan);
        countsmat(:,:,indChan,iter) = counts{thisChanName,iterNames(iter)};
    end
    counts{'avg',iterNames(iter)} = nanmean(countsmat(:,:,:,iter),3); % average over all channels for this iteration
end
counts.avg{'avg'} = nanmean(nanmean(countsmat,3),4); % grand average across all channels and iterations
counts.response_total('avg') = nanmean(counts.response_total);
counts.response_hz('avg') = nanmean(counts.response_hz);
counts.baseline_total('avg') = nanmean(counts.baseline_total);
counts.baseline_hz('avg') = nanmean(counts.baseline_hz);