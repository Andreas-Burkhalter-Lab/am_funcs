%% this version determine responsivity of a cell based on its most active trials rather than its preferred-stim-param trials

%%%% get_tuning_curves
% SP_files = spiking data for each plane generated by convertSuite2Pdata.m
% abf_file = behavioral + microscope timing file
% stimdata_file = data on stimulus parameters for each trial
% scopetiming_file = file listing start and stop points of scope within abf file generated with getAbfStartStop.m   
% regops_file = file generated during registration
%%%%% last updated 3/3/18
function res = get_tuning_curves_responsivityFromMostActiveTrials(SP_files, abf_file, stimdata_file, scopetiming_file, regops_file)

% for each cell, check this fraction of trials with the highest activity for responsivity
% to stimuli
frac_trials_to_check_for_responsivity = 1/7; 
% look in this time window in seconds for activity to compare stimulus-evoked and baseline activity 
responsivity_window_poststim_sec = [0 4] ; % [window start, window stop]
baseline_window_sec = 1.5; % use this number of seconds before stim onset for determining baseline

if ~exist('SP_files','var') || isempty(SP_files)
    SP_files = uigetfile('*.mat','Select ''_SP_'' .mat file(s).','MultiSelect','on');
end
if ~iscell(SP_files)
    SP_files = {SP_files};
end
nplanes_to_process = length(SP_files);

if ~exist('abf_file','var') || isempty(abf_file)
    abf_file = uigetfile('*.abf','Select .abf file.');
end

if ~exist('stimdata_file','var') || isempty(stimdata_file)
    stimdata_file = uigetfile('*.mat','Select stimdata .mat file.');
end

if ~exist('scopetiming_file','var') || isempty(scopetiming_file)
    scopetiming_file = uigetfile('*_scopetiming.mat','Select scopetiming .mat file.');
end

[abf_timepoints abf_sampling_interval_us abfinfo] = abfload(abf_file); % abf_sampling_interval_us = sampling interval microseconds
load(stimdata_file);
load(scopetiming_file);
load(regops_file); 
nplanes_in_session_total = ops1{1}.numPlanes;
ntrials = height(par_sets);
stimdur_us = pars.stimdur * 1e6; % convert from s to microseconds
stimdur_scans = round(stimdur_us / abf_sampling_interval_us);
responsivity_window_poststim_us = responsivity_window_poststim_sec .* 1e6; % convert from s to microseconds
responsivity_window_poststim_scans = round(responsivity_window_poststim_us ./ abf_sampling_interval_us); 
baseline_window_us = baseline_window_sec * 1e6;
baseline_window_scans = baseline_window_us / abf_sampling_interval_us; 

% get stim events from abf file
stimchan = find(strcmp(abfinfo.recChNames,'Psycho_In'));
stim_events = find_stimchan_events(abf_timepoints(:,stimchan)',pars,par_sets,abf_sampling_interval_us);
stim_events.Properties.VariableNames{find(strcmp('onset',stim_events.Properties.VariableNames))} = 'stim_onset';
stim_events.Properties.VariableNames{find(strcmp('duration',stim_events.Properties.VariableNames))} = 'stim_duration';
stim_events.chan_val = [];
stim_events.rspv_window_start = stim_events.stim_onset + responsivity_window_poststim_scans(1); % which scan to start the responsivity window on
stim_events.rspv_window_end = stim_events.stim_onset + responsivity_window_poststim_scans(2); % which scan to end the responsivity window on
stim_events.baseline_window_start = stim_events.stim_onset - baseline_window_scans; % which scan to start baseline window on

par_sets = [par_sets stim_events];

% get scope events from abf file
scope_events = find_stimchan_events(abf_timepoints(:,scopetiming.scope_chan)');

%% find stim events
%%% organize cells by plane - only processes 1 plane for now
for iplane = 1:nplanes_to_process
    load(SP_files{iplane})
    F_matrix = NaN(height(par_sets),size(F,2));
    par_sets = [par_sets table(F_matrix,F_matrix,F_matrix,'VariableNames',{'F_during_stim','F_during_rspv_window','F_prestim'})];
    sevents = scope_events(iplane:nplanes_in_session_total:end,:); % get scope events corresponding to this plane
    nimages = size(F,1);
    ncells = size(F,2);
    sevents = sevents(1:nimages,:); % ignore scope events that don't have a corresponding image
    sevents = [sevents table(F)];
    % for each stim trial, get corresponding scope events
    for itrial = 1:ntrials
        % F during stim
        sevents_in_range = find(sevents.onset>par_sets.stim_onset(itrial) & sevents.onset<[par_sets.stim_onset(itrial)+stimdur_scans]);
        sevents_in_range = sevents_in_range(1:end-1); % last image duration will be cut off by stim ending
        par_sets.F_during_stim(itrial,:) = mean(sevents.F(sevents_in_range,:));
        % F during stim 'responsivity window'
        sevents_in_range = find(sevents.onset > par_sets.rspv_window_start(itrial) & sevents.onset < par_sets.rspv_window_end(itrial));
        par_sets.F_during_rspv_window(itrial,:) = mean(sevents.F(sevents_in_range,:));
        % F during prestim 'baseline window'
        sevents_in_range = find(sevents.onset > par_sets.baseline_window_start(itrial) & sevents.onset < par_sets.stim_onset(itrial));
        sevents_in_range = sevents_in_range(1:end-1); % last image duration will be interrupted by stim onset 
        par_sets.F_prestim(itrial,:) = mean(sevents.F(sevents_in_range,:));
        
    end
end

% determine whether or not cell was responsive to stimuli
responsivity = table(false(ncells,1),NaN(ncells,1),cell(ncells,1),cell(ncells,1),cell(ncells,1),'VariableNames',{'rspv','rspv_pval','most_active_trial_inds','poststim_F','prestim_F'});
for icell = 1:ncells
    thiscell_respwindow_F_all = par_sets.F_during_rspv_window(:,icell); % F values for all trial within responsivity window
    thresh_F = quantile(thiscell_respwindow_F_all, 1-frac_trials_to_check_for_responsivity);
    top_trials = find(thiscell_respwindow_F_all > thresh_F); % find most active trials within responsivity window
    responsivity.prestim_F{icell} = par_sets.F_prestim(top_trials,icell);
    responsivity.poststim_F{icell} = par_sets.F_during_rspv_window(top_trials,icell);
    [responsivity.rspv(icell), responsivity.rspv_pval(icell)] = ttest(par_sets.F_prestim(top_trials,icell), par_sets.F_during_stim(top_trials,icell));
    responsivity.most_active_trial_inds{icell} = top_trials;
end
res.responsivity = responsivity;
res.par_sets = par_sets;
res.sevents = sevents;

%% tuning curves

% why are so many sf and tf tuned cells coming out not sgn responsive? check 'top trials' 
% preferred sf tf orient
% running modulation
% hwhm
% mean response timecourse with error bars
% check for OFF responses? see F_17203_2018-1-24_plane3_proc_SP_3.mat sftf trials 1-3, cell at {80,253}  


% get tuning curves
for thispar = {'sf','tf','orient'}
    do_analysis = 0; % default
    switch thispar{:} % check if the parameter was tested for
        case 'sf'
            if isfield(pars,'n_sfs') && pars.n_sfs > 1
                do_analysis = 1;
            end
        case 'tf'
            if isfield(pars,'n_tfs') && pars.n_tfs > 1
                do_analysis = 1;
            end
        case 'orient' %%% assume that if orient was tested, zero sfs were set and 1 tf was set while varying angle
            if isfield(pars,'nAngles') && pars.nAngles > 1 && ~isfield(pars,'n_sfs') && pars.n_tfs == 1
                do_analysis = 1;
            end
    end
    if do_analysis
        switch thispar{:}
            case 'sf'
                stimparvals = unique(par_sets.sfreq);
                stimparvals(stimparvals == pars.sf_fixed) = [];
                nparvals = pars.n_sfs;
                par_sets_varname = 'sfreq';
            case 'tf'
                stimparvals = unique(par_sets.tfreq);
                stimparvals(stimparvals == pars.tf_fixed) = [];
                nparvals = pars.n_tfs;
                par_sets_varname = 'tfreq';
            case 'orient'
                stimparvals = uniqe(par_sets.Angle);
                nparvals = pars.nAngles;
                par_sets_varname = 'Angle';
        end
        % make tuning table with stim par vals as rows, neurons as columns
        tuning = table(stimparvals,cell(nparvals,1),NaN(nparvals,ncells),'VariableNames',{'stimparval','resp','mean_resp'});
        for indpar = 1:nparvals
            thisparval = stimparvals(indpar);
            tuning.resp{indpar} = par_sets.F_during_stim(par_sets{:,par_sets_varname} == thisparval,:);
            tuning.mean_resp(indpar,:) = mean(tuning.resp{indpar});
        end
        
        % for each cell, make table of F response for each trial sorted by stim parameter and trial 
        tuning_by_cell = table(NaN(ncells,2),false(ncells,1),NaN(ncells,1),cell(ncells,1),cell(ncells,1),'VariableNames',{'centeryx','sgnf','anovap','trials','cellimage'});
        for icell = 1:ncells 
            cellimage = sparse(squeeze(logical(masks(icell,:,:))));
            tuning_by_cell.cellimage{icell} = cellimage;
            % get center of mass of cell image - use this location as an identifier for the cell
            [rows, cols] = ndgrid(1:size(cellimage,1), 1:size(cellimage,2));
            tuning_by_cell.centeryx(icell,1) = round(sum(sum(rows(cellimage))) ./ sum(sum(cellimage))); % y center from top
            tuning_by_cell.centeryx(icell,2) = round(sum(sum(cols(cellimage))) ./ sum(sum(cellimage))); % x center from left
            tuning_by_cell.trials{icell} = table(stimparvals,NaN(nparvals,1),NaN(nparvals,pars.repetitions),'VariableNames',{thispar{:},'resp_mean','resp'});
            for indpar = 1:nparvals
                match = find(par_sets{:,par_sets_varname} == stimparvals(indpar));
                tuning_by_cell.trials{icell}.resp(indpar,:) = [par_sets.F_during_stim(match,icell)]';
                tuning_by_cell.trials{icell}.resp_mean(indpar) = mean(tuning_by_cell.trials{icell}.resp(indpar,:));
            end
            tuning_by_cell.anovap(icell) = anova1(tuning_by_cell.trials{icell}.resp',[],'off'); % transpose resp so that responses are grouped by stim parameter value
            tuning_by_cell.sgnf(icell) = tuning_by_cell.anovap(icell) < 0.05;
        end
        tuning_by_cell.rspv = responsivity.rspv;
        tuning_by_cell.rspv_pval = responsivity.rspv_pval;
        
        % consolidate tuning tables
        switch thispar{:}
            case 'sf'
                res.sf_tuning = tuning;
                res.sf_tuning_by_cell = tuning_by_cell;
            case 'tf'
                if isfield(res,'sf_tuning_by_cell') % combine sf and tf tuning by cell if possible
                    keyvars = {'centeryx','cellimage','rspv','rspv_pval'}; % shared between sf and tf
                    [keyvarsfound,ind] = intersect(tuning_by_cell.Properties.VariableNames,keyvars);
                    res.sf_tuning_by_cell.Properties.VariableNames = cellfun(@(x)['sf_' x],res.sf_tuning_by_cell.Properties.VariableNames,'UniformOutput',0);
                    tuning_by_cell.Properties.VariableNames = cellfun(@(x)['tf_' x],tuning_by_cell.Properties.VariableNames,'UniformOutput',0);
                    tuning_by_cell.Properties.VariableNames(ind) = keyvarsfound; % rename key vars
                    res.sf_tuning_by_cell.Properties.VariableNames(ind) = keyvarsfound; % rename key vars
                    tuning_by_cell.cellimage = []; % non-char cell can't be a key var
                    res.sftf_tuning_by_cell = join(res.sf_tuning_by_cell,tuning_by_cell);
                    res = rmfield(res,'sf_tuning_by_cell');
                    res.tf_tuning = tuning;
                else
                    res.tf_tuning = tuning;
                    res.tf_tuning_by_cell = tuning_by_cell;
                end
            case 'orient'
                res.orient_tuning = tuning;
                res.orient_tuning_by_cell = tuning_by_cell;
        end
        clear tuning tuning_by_cell stimparvals nparvals par_sets_varname
    end
end
    
        