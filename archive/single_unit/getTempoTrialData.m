%%% analyze tempo trial from multiple files and write into results files        
% 1/12/17 edited on msi

function getTempoTrialData(filelist)
global tuning_summary_listfile

%% Parameters
% bad_trials_file = Excel file where first column lists .htb filenames and
%   second column lists the first bad trial in the adjacent .htb file;
%   analysis will be performed only on trials before this trial in that
%   file
%%% tuning_summary_listfile contains protocol names in the first column and the
%%% corresponding summary filename of the file containing results for all analyzed 
%%% cells on this protocol in the second column
% filelist = cell array of .htb filenames with and absolute/relative paths
%   to analyze
tuning_summary_listfile = 'C:\Users\AM\Documents\Lab\recordings\tuning_summary_list.xlsx';
psth_summary_file = 'C:\Users\AM\Documents\Lab\recordings\spike_timing_summary.mat';
topdir = 'C:\Users\AM\Documents\Lab\recordings';
bad_trials_file = 'C:\Users\AM\Documents\Lab\recordings\bad_trials';
display_results_figures = 0; % turn off to suppress opening new figure windows during analysis
display_command_line = 0; % turn off to suppress command line output from analysis scripts
loomingSpeedDependentSpikeWindow = 1; % if on, spike counting window for looming stim will depend on looming speed

% parameters for determining which files/analyses to run
onlyAnalyzeSpecifiedProtocols = 0; % if true, only analyze files with protocols matching the names below
    specifiedProtocols = {'LOOMING_SPEED'}; % {'GRATING_SPATIAL_FREQ'}; % {'SIZE_TUNING'}; % only has effect if onlyAnalyzeSpecifiedProtocols = true
try_Analysis_Switchyard = 1; % do analysis run through Analysis_Switchyard.m
    duplicate_Analysis_Switchyard_results = 1; % if false, check to see whether analysis results for a run exist in the summary file before analyzing
        duplicate_manual_check = 1; % check with user before each duplication (only matters if duplicate_Analysis_Switchyard_results == true)
try_psth_analysis = 1; % analyze for above-baseline responses through getSpikeTimes.m
  
%% setup
% if no filelist specified, call generate_filelist.m and run all files
% listed in cell_recordings_filename
if ~exist('filelist','var')
    [allfiles files_with_dirs] = generate_filelist; 
end
    
% set definitions
TEMPO_Defs;
ProtocolDefs; 
Path_Defs;
global pooling; % same setting as in TEMPO_GUI.M
pooling = 0; % same setting as in TEMPO_GUI.M
[~,bad_trials] = xlsread(bad_trials_file); % import bad trials list
for i = 1:size(bad_trials,1)
    if ~strcmp(bad_trials{i,1}(end-3:end),'.htb')
        bad_trials{i,1} = [bad_trials{i,1},'.htb']; % add extension if necessary
    end
end
origdir = pwd;
[~, summaryFileNames] = xlsread(tuning_summary_listfile);
psth_analysis_was_run = 0; % turn to true if psth analysis is run

% set inputs for Analysis_Switchyard
StartOffset = 0;
StopOffset = 0; 
SpikeChan = 1;
SpikeChan2 = 1;
StartCode = 4; % default for most protocols
StopCode = 5; % default for most protocols
batch_flag = 0;  % tells analysis switchyard whether or not to process batch files
UseSyncPulses = 0;

%% load and analyze data from each file
nfiles = length(filelist);
for i = 1:nfiles
    if ~strcmp(filelist{i}(end-3:end),'.htb')
        filelist{i} = [filelist{i},'.htb']; % add extension if necessary
    end
end
for fileind = 1:nfiles
    cd(topdir);
    thisfile = filelist{fileind};
%     fprintf([thisfile '\n'])
    [PATH fname ext] = fileparts(thisfile);
    FILE = [fname ext];
    % select and load tempo data into matlab
    cd(PATH);
    PATH = [];
    
    % get protocol name
    [return_value, info_text] = GetHTBInfo(PATH,FILE);
    protocol_name = info_text{1}(18:end-1);
    if strcmp(protocol_name,'COHERENCE') % use the correct keyword
        protocol_name = 'MOTION_COHERENCE';
    end
    protocol_name_reformat = strrep(protocol_name,'_',' '); %replace underscore with space
    Protocol = find(strcmpi(protocol_names,protocol_name_reformat),1)-1; % subtract 1 because indices start with zero, not 1
    
    %% checks to decide whether or not to analyze this run
    if try_Analysis_Switchyard
    do_Analysis_Switchyard = 1; % for overwrite checking; start on by default
        if onlyAnalyzeSpecifiedProtocols
            if any(strcmp(protocol_name,specifiedProtocols)) % if this run is not one of the specified protocols
                do_Analysis_Switchyard = 1; % proceed to further checks
            elseif ~any(strcmp(protocol_name,specifiedProtocols)) % if the run is of a specified protocol to analyze
                do_Analysis_Switchyard = 0; % skip further checks, do not analyze this run
            end
        end

        % overwrite checks
        if do_Analysis_Switchyard && duplicate_Analysis_Switchyard_results && ~duplicate_manual_check % if automatically overwriting, save time by not opening summary file
            do_Analysis_Switchyard = 1; 
        elseif do_Analysis_Switchyard % if we are checking for overwrite and do_analysis hasn't already been blocked
            matchrow = find(strcmp(protocol_name,summaryFileNames(:,1)));
            if isempty(matchrow) % if listed file is not an implemented protocol yet
                warning('Protocol %s from file %s not yet implemented; will not analyze.', protocol_name, FILE);
                do_Analysis_Switchyard = 0;
            elseif ~isempty(matchrow) && exist(summaryFileNames{matchrow,2},'file'); % if a summary file for this protocol exists and it's an implemented protocol
                thisSummaryFile = summaryFileNames{matchrow,2};
                thisSummaryList = importdata(thisSummaryFile);
                % determine whether results for this run already exist in the summary file
                if any(find(strcmp(FILE,thisSummaryList.textdata))) % if analysis for this run already exists in the summary file
                    if duplicate_Analysis_Switchyard_results && duplicate_manual_check
                        overwrite_check_answer = input(sprintf('''%s'' is already listed in summary file ''%s''. Enter ''y'' to duplicate.\n',getfname(thisfile),getfname(thisSummaryFile)),'s');
                        if strcmp(overwrite_check_answer,'y')
                            do_Analysis_Switchyard = 1;
                        else
                            do_Analysis_Switchyard = 0;
                            fprintf('\nWill not overwrite file %s.',thisfile);
                        end
                    elseif ~duplicate_Analysis_Switchyard_results
                        do_Analysis_Switchyard = 0; % automatically do not overwrite
                    end
                end
            elseif ~exist(thisSummaryFile,'file') % summary files does not exist, create a new one and do analysis
                do_Analysis_Switchyard = 1; 
            end
        end
    end
            
    %% load data
    if try_psth_analysis || (try_Analysis_Switchyard && do_Analysis_Switchyard)
        % determine trial range
        [~, good_data] = LoadTEMPOData(PATH,FILE); % load data
        select_data = good_data; % analyze good data
        ntrials = size(good_data.spike_data,3); % total trials recorded
        bad_trial_list_index = find(strcmp(filelist{fileind},bad_trials(:,1))); % check if this trial is listed in bad_trials_file
        if ~isempty(bad_trial_list_index) % if this file has a bad trial, analyze only trials before the bad trial
            skip_trials = str2num(bad_trials{bad_trial_list_index,2});
        elseif isempty(bad_trial_list_index)
            skip_trials = [];
        end
        BegTrial = 1;
        EndTrial = ntrials; % handle bad trials through 'skip_trials', not BegTrial and EndTrial
    end
    
    %% run Analysis_Switchyard if applicable
    if try_Analysis_Switchyard && do_Analysis_Switchyard
        % choose which analysis to do; see 'analysis_strings' from ProtocolDefs.m
        implemented_protocol = 1; % will be turned off if we detect  protocol not yet implemented

        switch protocol_name 
            case 'RF_MAPPING' %%% rf mapping with drifting dots, not gratings
                Analysis = 'Fit with 2D Gaussian';
            case 'GRATING_RF_MAP'
                Analysis = 'Fit 2D Gaussian';
        %         Analysis = 'Plot PSTH';
        %         Analysis = 'Plot Spike Rasters';
            case 'SIZE_TUNING' %%% dot or dot field size, not grating size
                Analysis = 'Plot Tuning Curve';
    %             Analysis = 'Plot Rasters/Histograms';
    %             Analysis = 'Plot Vergence Data';
    %             Analysis = 'Plot Auto and Cross Correlograms';
            case 'GRATING_SIZE'
                Analysis = 'Plot Tuning Curve';
    %             Analysis = 'Plot PSTH';
    %             Analysis = 'Plot Spike Rasters';
    %             Analysis = 'Pooling';
            case 'GRATING_SPATIAL_FREQ'
                Analysis = 'Plot Tuning Curve';
%                 Analysis = 'Plot PSTH';
        %         Analysis = 'Plot Spike Rasters';
        %         Analysis = 'PSTH Fourier Analysis';
        %         Analysis = 'Stimulus Onset';
            case 'GRATING_TEMPORAL_FREQ'
                Analysis = 'Plot Tuning Curve';
        %         Analysis = 'Plot PSTH'
        %         Analysis = 'Plot Spike Rasters'
            case 'SPEED_TUNING' % drifting dot field speed
                Analysis = 'Plot Tuning Curve';
    %             Analysis = 'Plot Histograms';
    %             Analysis = 'Plot Event Times';
    %             Analysis = 'Plot PSTH';
    %             Analysis = 'Plot Auto Correlograms';
    %             Analysis = 'Plot Cross Correlograms';
    %             Analysis = 'Plot ISI Histogram';
    %             Analysis = 'Plot Joint PSTH';
    %             Analysis = 'Experiment Playback';
    %             Analysis = 'Fit LogGauss (Harris)';
    %            Analysis =  'Fit Gaussian';
            case 'MOTION_COHERENCE' % drifting dot field motion coherence
                Analysis = 'Plot Tuning Curve';
    %             Analysis = 'Plot PSTH';
    %             Analysis = 'Plot Spike Rasters';
            case 'LOOMING_SPEED'
                Analysis = 'Plot Tuning Curve';
        %         Analysis = 'Plot Histograms';
        %         Analysis = 'Plot Auto and Cross Correlograms';
        %         Analysis = 'Fit Gaussian';
            case 'GRATING_ORIENTATION'
                Analysis = 'Plot Tuning Curve';
        %         Analysis = 'Plot PSTH'
        %         Analysis = 'Plot Spike Rasters'
        %         Analysis = 'Stimulus Onset'
            case 'DIRECTION_TUNING' 
                Analysis = 'Plot Tuning Curve';
        %         Analysis = 'Plot Event Times';
        %         Analysis = 'Plot Cross Correlograms';
        %         Analysis = 'Plot Auto Correlograms';
        %         Analysis = 'Plot ISI Histogram';
        %         Analysis = 'Plot PSTH';
        %         Analysis = 'Plot Spike Rasters';
        %         Analysis = 'Experiment Playback';  
            case 'OPTIC_FLOW_SPEED'
                Analysis = 'Plot Tuning Curve';
        %         Analysis = 'Plot Histograms';
        %         Analysis = 'Plot Auto and Cross Correlograms';
        %         Analysis = 'Fit Gaussian';
            case 'OPTIC_FLOW_COHER'
                Analysis = 'Plot Tuning Curve';
        %         Analysis = 'Plot PSTH'
        %         Analysis = 'Plot Spike Rasters'
            otherwise 
                implemented_protocol = 0;
                error('Protocol ''%s'' used in %s is not yet implemented in this analysis script.',...
                    protocol_name,thisfile)
        end

        if implemented_protocol
            original_visibleFigure_value = get(0,'DefaultFigureVisible');
            if display_results_figures
                set(0,'DefaultFigureVisible','on');
            else
                set(0,'DefaultFigureVisible','off');
            end

            if ~isempty(skip_trials) && ~any(strcmp(Analysis,{'Fit with 2D Gaussian','Fit 2D Gaussian','Plot Tuning Curve'}))
                error(sprintf(['This run contains bad trials, which are only implemented to be processed with the following analysis types:'...
                    '\n Fit with 2D Gaussian \n Fit 2D Gaussian \n Plot Tuning Curve']))
            end
             
            %% run protocol-specific analysis
            if display_command_line
                Analysis_Switchyard(select_data, Protocol, {Analysis}, SpikeChan, SpikeChan2, StartCode, StopCode,...
                BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses, skip_trials);
            else % suppress command line output
                evalstring = ['Analysis_Switchyard(select_data, Protocol, {Analysis}, SpikeChan, SpikeChan2, StartCode, StopCode,',...
                'BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses, skip_trials);'];
                evalc(evalstring);
            end
            

            set(0,'DefaultFigureVisible',original_visibleFigure_value) % reset figure display mode
        end
    end
    %% run getSpikeTimes if applicable
    if onlyAnalyzeSpecifiedProtocols
        if any(strcmp(protocol_name,specifiedProtocols)) % if this run is not one of the specified protocols
            do_psth_analysis = 1; % proceed to further checks
        elseif ~any(strcmp(protocol_name,specifiedProtocols)) % if the run is of a specified protocol to analyze
            do_psth_analysis = 0; % skip further checks, do not analyze this run
        end
    elseif ~onlyAnalyzeSpecifiedProtocols
        do_psth_analysis = 1;
    end
    if try_psth_analysis && do_psth_analysis
        num_trials = size(select_data.event_data, 3);
        
        [StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin] = CheckTimeOffset(select_data, num_trials, StartCode, StopCode, StartOffset, StopOffset, UseSyncPulses);
        if loomingSpeedDependentSpikeWindow && any(Protocol == [LOOMING_SPEED; LOOMING_GAMP; LOOMING_GASSIAN_AMP]) % if it's a looming protocol and loomingSpeedDependentSpikeWindow==1
            select_data.spike_rates = ComputeSpikeRates_looming(select_data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin); % use spike window based on looming speed        
        elseif ~loomingSpeedDependentSpikeWindow || ~any(Protocol == [LOOMING_SPEED; LOOMING_GAMP; LOOMING_GASSIAN_AMP]); % if it's not looming or we are using fixed spike windows for looming
            select_data.spike_rates = ComputeSpikeRates(select_data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin); % use same spike window for all stim params
        end
        
        getSpikeTimes(thisfile, FILE, psth_summary_file, protocol_name, select_data, SpikeChan, StartCode, StopCode, BegTrial,...
            EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, pooling, skip_trials, loomingSpeedDependentSpikeWindow); 
        psth_analysis_was_run = 1;
    end
end

% tabulate spike timing results, save as table
if psth_analysis_was_run
    analyze_spike_timing(loomingSpeedDependentSpikeWindow);
end

clear global tuning_summary_listfile
cd(origdir)