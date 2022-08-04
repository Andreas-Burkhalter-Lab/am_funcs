%GETSPIKES: Get spikes from raw .merec file to align with stimchan events.
%%% Last updated 1/7/16 on msi


%% Snip and sort if necessary
do_autosnip2_autosort = 1; % initialize
if ~rerun_autosnip2_autosort
    if exist([datadir filesep 'overview.mat'],'file')
        load([datadir filesep 'overview.mat'],'sorthead');
        if strcmp(getfname(sorthead.fh.filename), getfname(opt_autosnip2.files{:}))
            disp('Found autosort overview file; not rerunning autosnip2 and autosort.')
            do_autosnip2_autosort = 0;
        end
    end
end
if do_autosnip2_autosort
    % Filtering
    if do_filter 
        passband = [highpass 0.48*merec_obj.scanrate]; % highpass slightly below Nyquist
        opt_snippetfile = snipoptions(struct(... % make the conditioning filter
            'Fs',merec_obj.scanrate,'bandpass',passband)); 
    else % use snippetfile defaults only... something goes wrong if we filter in snippetfile and not in autosnip2
        opt_snippetfile = []; 
        clear passband highpass
    end

    autosnip2_AM(opt_autosnip2, opt_snippetfile);
    autosort_AM(snipfile2sortheader(strcat(opt_autosnip2.files{:}(1:end-6),'.ssnp')), datadir);
    if ~do_filter
        clear opt_snippetfile % clear if irrelevant
    end
else % clear irrelevant variables
    clear opt_autosnip2 snip_range_ms do_filter passband opt_snippetfile highpass
end
load([datadir '/overview.mat'])

%% Get snips from .ssnp and autosort files
% Find the correct snip file.
if exist(fullfile(sorthead.fh.abspathstr,sorthead.fh.filename),'file') % if the .ssnp file is in original location
    snip_file = fullfile(sorthead.fh.abspathstr,sorthead.fh.filename);
elseif exist([getfname(resp_file) '.ssnp'], 'file')    
    snip_file = [getfname(resp_file) '.ssnp'];
else 
    snip_file = uigetfile('*.ssnp','Select snip .ssnp file.');
end
    
snipdata = sorthead_from_raw(sorthead, opt_autosnip2.snip_channels, snip_file);
snipdata.unitnames = strrep(cellfun(@(x)['ch',num2str(x)],num2cell(rec_chans),'UniformOutput',false),'.','_');

% maybe list on a single line which unit is being worked on by sorted_from_raw  
%%% or find a way to speed up sorted_from_raw
% modified snipdata.channels field not yet tested; supposed to list the channel 
%   corresponding to each cluster, to make it compatible with manual_trial_noise_check 
if manual_clustering 
    cass_sort_apply_AM(datadir,[],snip_file); 
    tempsnips = {}; tempsniptimes = {}; tempunitnames = {};
    chansToAppend = []; % append to front of snipdata.channels
    for indChan = 1:length(snipdata.unitnames)
        if exist([datadir filesep 'chan' num2str(rec_chans(indChan))],'dir')  
            tempdir = dir([datadir filesep 'chan' num2str(rec_chans(indChan))]);
            [junk junk ext] = cellfun(@fileparts,extractfield(tempdir,'name'),'UniformOutput',false);
            issortfile = strcmp(ext,'.sorted');
            chanName = rec_chans(indChan); % channel label from .merec file
            if length(find(issortfile)) > 1
                error('There must be no more than 1 ''.sorted'' file in each chan diretory.')
            elseif length(find(issortfile)) == 1; % if there is one cass sorting file for this chan
                sortedfile = [datadir filesep 'chan' num2str(chanName) filesep tempdir(issortfile).name];
                thisChClusts = sorted_from_raw(merec_obj, sortedfile, opt_autosnip2.snip_range); % may take a long time
                nClustsThisChan = length(thisChClusts);
                tempsniptimes = [(flipud(importdata(sortedfile,'-mat')))' tempsniptimes];
                for clust = 1:nClustsThisChan
                    tempsnips = [thisChClusts(clust) tempsnips]; % preallocate for speed?
                    tempunitnames = [['ch' num2str(chanName) '_cl' num2str(clust)] tempunitnames];
                end
                
                %%% add repeats to snipdata.channels to correspond to
                %%% snipdata.unitnames - NOT YET TESTED
                chansToAppend = [chansToAppend repmat(chanName,1,nClustsThisChan)];
            end
        end
    end
    % Add manually clustered snip data to unclustered snip data. 
    snipdata.snips = [fliplr(tempsnips) snipdata.snips]; 
    snipdata.sniptimes = [fliplr(tempsniptimes) snipdata.sniptimes];
    snipdata.unitnames = [fliplr(tempunitnames) snipdata.unitnames];
% % % % %     snipdata = rmfield(snipdata,'channels'); % remove to avoid confusion after we add more unit -REPLACED WITH LINE BELOW
    snipdata.channels = [chansToAppend snipdata.channels]; % NOT TESTED
end

% If only analyzing a limited range of trials, cut down trialdata.
%%% would be faster to find the scan range of the last trial to analyze,
%%% then limit importing of cluster spike times to this scan range. 
if analyzeLimitedTrialRange
    fprintf('Only analyzing trials %g through %g.',...
        trialStartEndToAnalyze(1), trialStartEndToAnalyze(2))
    trialdata = trialdata(trialStartEndToAnalyze(1):trialStartEndToAnalyze(2),:);
else
    clear trialStartEndToAnalyze
end

% Get spike counts for each event for each unit.
% trialdata.latency contains latency in seconds relative to stim onset;
%   baseline latencies are negative, response latencies are positive. 
% trialdata.snip_ind contains the indices of the spikes found on this unit
%   in this event block. To get their waveforms, use:
%   snipdata.snips{unit}(:,trialdata.snip_ind{event,unit})
%%% Rather than getting spikes from the actual stim duration in
%%%     trialdata.stimdur, get spikes during an fixed-duration time period for each
%%%     trial. Baseline spike data is also taken from fixed-duration
%%%     pre-stim period, rather than actual time between previous offset
%%%     and current onset. 
ntrials = size(trialdata,1);
nunits = length(snipdata.unitnames);
isi_scans = ss_pars.isi * merec_obj.scanrate; 
trialdata.spikes = NaN(ntrials,nunits);
trialdata.baseline_spikes = NaN(ntrials,nunits);
trialdata.spike_latency = cell(ntrials,nunits);
trialdata.snip_ind = cell(ntrials,nunits);
trialdata.baseline_snip_ind = cell(ntrials,nunits);
trialdata.baseline_spike_latency = cell(ntrials,nunits);

for trial = 1:ntrials
    if disp_trialBeingAnalyzed
        fprintf('\nAnalyzing trial %g/%g...',trial,ntrials)
    end
    thisTrialOnset = trialdata.onset(trial);
    for thisunit = 1:nunits
        thisUnitSnipTimes = snipdata.sniptimes{thisunit};
        trialdata.snip_ind{trial,thisunit} = find(thisUnitSnipTimes > thisTrialOnset...
            & thisUnitSnipTimes <= thisTrialOnset+window_scans); % response snip inds
        trialdata.baseline_snip_ind{trial,thisunit} = find(thisUnitSnipTimes >= thisTrialOnset-isi_scans...
            & thisUnitSnipTimes < thisTrialOnset); % prestim baseline snip inds
        trialdata.spike_latency{trial,thisunit} = ...
            double(thisUnitSnipTimes(trialdata.snip_ind{trial,thisunit}) - thisTrialOnset) ./ merec_obj.scanrate; % secs
        trialdata.baseline_spike_latency{trial,thisunit} = ...
            double(thisUnitSnipTimes(trialdata.baseline_snip_ind{trial,thisunit}) - thisTrialOnset) ./ merec_obj.scanrate; % secs
        trialdata.spikes(trial,thisunit) = length(trialdata.snip_ind{trial,thisunit});
        trialdata.baseline_spikes(trial,thisunit) = length(trialdata.baseline_snip_ind{trial,thisunit});
    end
end

% Manual trial checking
if manual_trial_check % add variable for spikes automatically detected, with manually removed trials turned to NaNs
    [flagged_trials trialnoise nFlaggedTrials] = manual_trial_noise_check(trialdata,snipdata,opt_trialCheck,merec_obj,window_scans);
    trialdata.discard = trialnoise.flagged; % if =1 for a given trial, this trial will not be counted
    trialdata = trialdata(:,{'discard','diam','Angle','sf','tf','onset','dur','spikes','snip_ind',...
        'spike_latency','baseline_spikes','baseline_snip_ind','baseline_spike_latency'}); % reorder variables
    if opt_trialCheck.keepRaw % save spike/trial data before removing noisy trials
        trialdata_raw = trialdata;
    end
    trialdata.spikes(trialdata.discard==1,:) = NaN; % set spike data to NaN for all trials flagged as noisy
end