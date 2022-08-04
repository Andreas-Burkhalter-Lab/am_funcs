% getSpikesTimes: saves into a .mat file the spike times from trials in which the optimal stimulus was
% presented, including from a pre-stim period
%%%% adapted from param_by_trial_PSTH.m by AM
%%%% last updated 8/29/16 
%%% datafile is filename and extension of run to be analyzed
%%% full_file is dir+filename+extension of run
%%% summaryfile is a .mat file in which spike timing data from all runs
%%%    will be saved

%Edited on 1/20/08 by RLS to incorporate a seperate 'Pooling' function for use
%with calculating latency
%Edited on 10/05/07 by RLS to display file name on the top of the plots
%Edited on 06/19/06 by CMA to remove all references to SpikeChan2 or
%spikes2, which produce the error "??? Undefined function or variable
%'SpikeChan2'." when run with our data.  The original file is backed up in
%the folder Chris's Lab Tools Backup located on the desktop of this
%machine.

%%%% 8/13/16 AM added conditional when writing outfile to write into the correct
%%%%    directory if on a specific computer; functionality should be
%%%%    unchanged if using on any lab computer

function getSpikeTimes(datafile, full_file, summaryfile, protocol_name, data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial,...
    StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Pooling, skip_trials, loomingSpeedDependentSpikeWindow)

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;
symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'm-' 'g-' 'b-' 'b--' 'r--' 'g--' 'c--'};

% matches protocol names to stim param keyword numbers provided by
% ProtocolDefs; need to find stim parameter keyword for determining optimal
% stim param value
%%% third column instructs to look for parameter in data.dots_params vs.
%%% data.gratings_params
protocol_keyword_pairs = {...
    'SIZE_TUNING', DOTS_AP_XSIZ, 'dot';... % for dot size tuning 
    'GRATING_SIZE', GRAT_WIDTH, 'grat';... 
    'GRATING_SPATIAL_FREQ', GRAT_SPATIAL_FREQ, 'grat';...
    'GRATING_TEMPORAL_FREQ', GRAT_TEMPORAL_FREQ, 'grat';... 
    'SPEED_TUNING', DOTS_SPEED, 'dot';...
    'MOTION_COHERENCE', DOTS_COHER, 'dot';...
    'LOOMING_SPEED', DOTS_AP_VEL, 'dot';... 
    'DIRECTION_TUNING', DOTS_DIREC, 'dot';...
    'GRATING_ORIENTATION', GRAT_ORIENTATION, 'grat';...
    'OPTIC_FLOW_COHER', DOTS_COHER, 'dot';...
    'OPTIC_FLOW_SPEED', DOTS_AP_VEL, 'dot'...
    };

col = find(strcmp(protocol_name,protocol_keyword_pairs(:,1))); 
if isempty(col)
    warning('Protocol %s in file %s not yet implemented for batch psth analysis',protocol_name,datafile)
    implemented_prot = 0;
elseif ~isempty(col)
    implemented_prot = 1;
    paramNumber = protocol_keyword_pairs{col,2};
    dot_or_grat = protocol_keyword_pairs{col,3};
end

% check to see whether or not this file already has an entry in the summary file     
do_analysis = 0;
if implemented_prot
    if ~exist(summaryfile,'file')
        do_analysis = 1;
    elseif exist(summaryfile,'file')
        mfile = matfile(summaryfile,'Writable',true);
        filelist = mfile.spikes_list(:,1); % this column of the cell array should be the list of files already analyzed
        if any(strcmp(filelist,full_file))
            add_duplicate = input(sprintf('''%s'' is already listed in summary file ''%s''. Enter ''y'' to add duplicate entry.\n',datafile,getfname(summaryfile)),'s');
            if strcmp(add_duplicate,'y')
                do_analysis = 1;
            elseif ~strcmp(add_duplicate,'y')
                do_analysis = 0; % don't add duplicate entry 
            end
        else % entry does not yet exist for this run
            do_analysis = 1;
        end
    end
elseif ~implemented_prot
    do_analysis = 0;
end
    
if do_analysis
    if strcmp(dot_or_grat,'dot')
        param_by_trial = data.dots_params(paramNumber,:,PATCH1);
    elseif strcmp(dot_or_grat,'grat')
        param_by_trial = data.gratings_params(paramNumber, :, PATCH1);
    end
        %get the column of values of param_by_trials in gratings_params matrix
        %     error if all param vals are the same
    if length(unique(param_by_trial)) < 2
        error(sprintf('Only found 1 parameter value for %s across all trials in file %s.',protocol_name,full_file))
    end
    
    %get indices of any NULL conditions (for measuring spontaneous activity)
    null_trials = logical( (param_by_trial == data.one_time_params(NULL_VALUE)) );
    clear has_null_trials
    if ~isempty(find(null_trials)) % indicate whether or not null trials are available for baseline data
        has_null_trials = true;
    elseif isempty(find(null_trials))
        has_null_trials= false;
    end
    
    %get the column of stimulus type values

    stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
    unique_stim_type = munique(stim_type(~null_trials)');

    %now, get the firing rates for all the trials
    spike_rates = data.spike_rates(SpikeChan, :);

    unique_param_by_trial = munique(param_by_trial(~null_trials)');

    %now, remove trials that do not fall between BegTrial and EndTrial
    trials = 1:length(param_by_trial);		% a vector of trial indices
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
    
    % mark bad trials to skip
    if exist('skip_trials','var')
        select_trials(skip_trials) = false;
    end
    
    count=0;
    nullcount=0;
    %subplot(2, 1, 2);

    for j = 1:length(unique_stim_type);

        grating_select = logical( (stim_type == unique_stim_type(j)) );
        plot_x=param_by_trial(~null_trials & select_trials & grating_select);
        plot_y=spike_rates(~null_trials & select_trials & grating_select);

        %NOTE:  Inputs to PlotTuningCurve.m must be column vectors, not row vectors
        %because of munique().
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
        %NOTE: last argument=0 means just get output, no plot

        %py is the response and px is the stimulus value.  We are interested in
        %the cell's response at optimized conditions, therefore, we find the
        %max of py, and select the corresponding px as stimulus value that
        %results in optimized response.
        [~, index] = max(py);
        optimized_stimulus = px(index);

        optTrials = find(param_by_trial == optimized_stimulus); % trials with the optimal stimulus
        spikes_opt = cell(length(optTrials),1); % each cell contains the spike times for one trial
        spikes_null = cell(length(find(null_trials)),1); % will be empty cell if there are no null trials
        for trial=1:length(param_by_trial);
            % optimal stimulus
            if (param_by_trial(trial) == optimized_stimulus & stim_type(trial) == unique_stim_type(j));
                count = count+1;
                spikes =  find(data.spike_data(SpikeChan,:,trial) == 1);
                if j == 1;
                    spikes_opt{count} = spikes; 
                elseif j == 2;
                    error('Error: spike chan 2 not yet supported')
                else
                    error('Error: More than two grating types')
                end
                hold on;
            % null trials
            elseif has_null_trials && param_by_trial(trial)==-9999 && stim_type(trial) == unique_stim_type(j);
                nullcount=nullcount+1;
                nullspikes = find(data.spike_data(SpikeChan,:,trial) == 1);
                spikes_null{nullcount} = nullspikes;
            end
        end
        trials(j) = count-1;
    end
    
    % if looming run and if loomingSpeedDependentSpikeWindow == 1, get
    % response rates for optimized recorded by ComputeSpikeRates_looming.m
    % in data.spike_rates; assume we are only using channel 1
    if loomingSpeedDependentSpikeWindow && any(Protocol == [LOOMING_SPEED; LOOMING_GAMP; LOOMING_GASSIAN_AMP]) % if it's a looming protocol and loomingSpeedDependentSpikeWindow==1
        spikes_loom_rsp = data.spike_rates(1,optTrials); 
    elseif ~loomingSpeedDependentSpikeWindow || ~any(Protocol == [LOOMING_SPEED; LOOMING_GAMP; LOOMING_GASSIAN_AMP]); % if it's not looming or we are using fixed spike windows for looming
        spikes_loom_rsp = [];
    end

    % save results into running results file of spike timing results
    newEntry = {full_file,spikes_opt',spikes_null',has_null_trials,StartEventBin(trial),StopEventBin(trial),protocol_name,spikes_loom_rsp}; % spike data for this run
    if exist(summaryfile,'file')
        arrayHeight = size(mfile.spikes_list,1);
        mfile.spikes_list(arrayHeight+1,:) = newEntry; %% add spike timing data for this file to the running table of spike timing results
    elseif ~exist(summaryfile,'file')
        spikes_list = newEntry;
        column_headings = {'file_name','spikes_opt','spikes_null','has_null_trials','stim_onset_time','stim_offset_time','protocol_name','spikes_loom_rsp'}; % column labels for the spikes_list array
        save(summaryfile,'spikes_list','column_headings','-v7.3'); % create new file in which to save spike timing results
    end
end