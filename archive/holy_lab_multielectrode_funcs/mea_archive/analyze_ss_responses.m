function [] = analyze_ss_responses(datadir, resp_file, stim_file, rec_chans)
%ANALYZE_SS_RESPONSES Create tuning curves surround-suppression response
%data.
%   resp_file is the .merec file of MEA responses
%   stim_file is the stimulus log file created on the stimulus computer
%   rec_chans specifies which channels contain electrode data to analyze
%       (default = all channels except 'reserved_stimchans')
%   datadir specifies the directory in which to save results and look for
%       autosnip2, autosort, and cass results (default = current directory) 
%%% See surround_suppression_stim and run_experiment for relevant units. 
%%% should we subtract baseline from response firing rates?
%%% Last updated 7/13/15 on vivid

% how to get gaussian distr of peak amps within a cluster? 


stim_chan = 63;
resp_window = 2; % post=stim time (s) in which to count spikes as responding to the stim
rerun_autosnip2_autosort = 0;
manual_clustering = 1; % run cass_sort_apply and treat all manual clusters as separate channels
reserved_stimchans = [55 63]; % do not treat these channels as electrode channels by default
snip_range_ms = [-1 2]; % spike-relative snippet window in ms (default = [-1 3])

opt_autosnip2.thresh_factor = 1; % =2 in Vaiceliunaite et al. 2013: thresh  = twice 6*median
opt_autosnip2.resume = 0;
opt_autosnip2.feedback_channels = 63;
opt_autosnip2.do_filter = 0;
opt_autosnip2.compare_hz60 = 0;
opt_autosnip2.time4data_segment = 5;
opt_autosnip2.polarity = -1;
opt_autosnip2.do_in_time_order = 1;
opt_autosnip2.extension = '.ssnp';
opt_autosnip2.iterative_threshold = 0; 
opt_autosnip2.no_env = 1; 
opt_autosnip2.no_redo_vlv = 1;

if ~exist('resp_file','var') || isempty(resp_file)
    resp_file = uigetfile('*.merec','Select response file.');
end
opt_autosnip2.files = {resp_file};
merec_obj = merecmm(resp_file);
opt_autosnip2.snip_range = merec_obj.scanrate*1e-3.*snip_range_ms;
window_scans = resp_window*merec_obj.scanrate;

if ~exist('stim_file','var') || isempty(stim_file)
    stim_file = uigetfile('*.mat','Select stimulus log file.');
end
load(stim_file);

if ~exist('rec_chans','var') || isempty('rec_chans')
    disp(sprintf(['Using all channels except [%s].'], num2str(reserved_stimchans)));
    allchans = merec_obj.channels;
    rec_chans = allchans(arrayfun(@(x)(~any(x==reserved_stimchans)), allchans));
end
opt_autosnip2.snip_channels = rec_chans;    

if ~exist('datadir','var') || isempty(datadir)
    datadir = pwd;
end

%% Extract event timing
int_wave = round(merec_obj([stim_chan], [1:merec_obj.nscans]));  %% round stim wave to nearest integer
transitions(1,:) = [1, find(~diff(int_wave)==0) + 1]; % event timings, starting with initial state (1)... +1 is to find post-transition time
transitions(2,:) = int_wave(transitions(1,:)); %% get stim channel values after transitioning from one value to another   
transitions(3,:) = [diff(transitions(1,:)) NaN]; % event durations  (ignore last nonstim event duration)

% Delete all event transitions immediately following the preceding transition, 
% in case it takes more than one timestep to complete a transition. 
trans_diff = [diff(transitions(1,:)) 0]; %% get time length in time steps between events, append 0 for last timestep
count = 0;
for i = 1:size(transitions,2)
    if trans_diff(i)==1 %% if this transition is 1 timestep before the next transition
        count = count+1; 
    elseif count  %% if this transition is 1 step away from the preceding but not the following transition
        transitions(1,i-count+1:i) = -1;        %% mark problematic transitions with -1 for deletion
        transitions(2,i-count) = transitions(2,i); % use the final stimchan value in this set of adjacent transitions
        transitions(3,i-count) = transitions(3,i)+count; %% correct the duration of nondeleted event
        count = 0;
    end
end
transitions = transitions(:,transitions(1,:)>0);  %% eliminate all transitions one timestep after the preceding transition
    
% If the initial pre-stim Vstimoff value is preceded by another stimchan
% value, give a warning but delete the first value. The first value was
% probably the stimchan value before stimulus presentation began. 
if transitions(2,2)==pulse.Vstimoff & transitions(2,1)~=pulse.Vstimoff
   warning(['Ignoring stimchan value preceding first stim-off channel value'...
       ' (value = %g from scan 1 to scan %g).'],transitions(end,1),transitions(1,2));
   transitions = transitions(:,2:end); % delete the preceding nonzero column
end

% Check that stim events follow the prescribed order of 
% Vstimoff, Vstimon, Vstimoff...Vstimoff
if rem(length(transitions)-1,2) % if 'transitions'-1 is NOT divisible by 2
    error('Event timing-parsing error: number of events must equal 2 * nPresentations + 1')
elseif any(transitions(2,1:2:end) ~= pulse.Vstimoff) ||...    % if Vstimon values are inccorrect
       any(transitions(2,2:2:end) ~= pulse.Vstimon)      % ... or Vstimoff values are incorrect
   error('Event timing-parsing error: events must follow this sequence: [Vstimoff,Vstimon,Vstimoff... VStimoff].')
end

transitions(:,1:2:end) = [];    %% eliminate stim-absent (Vstimoff) events

% Check that response window is shorter than the time between stim onsets. (This will not the check final stim vs. end of recording.) 
if any(diff(transitions(1,:) < window_scans))
    error('Event timing-parsing error: response window to analyze must be shorter than time between stimulus onsets.')
end

if size(transitions,2) == size(par_sets,1)
    par_sets.onset = (transitions(1,:))'; % stim onset time in scans
    par_sets.dur = (transitions(3,:))'; % duration in scans
    clear transitions trans_diff count
else % if number of stim chan events does not match number indicated on stim comp
    error(['Number of events found on stim channel (%g) does not equal number of '...
     'events recorded on stimulus computer (%g)'],size(transitions,2),size(par_sets,1))
end

%% Get spike count data
do_autosnip2_autosort = 1;
if ~rerun_autosnip2_autosort
    if exist([datadir filesep 'overview.mat'],'file')
        load([datadir filesep 'overview.mat'],'sorthead');
        if strcmp(sorthead.fh.filename(1:end-5), opt_autosnip2.files{:}(1:end-6))
            disp('Found autosort overview file; not rerunning autosnip2 and autosort.')
            do_autosnip2_autosort = 0;
        end
    end
end
if do_autosnip2_autosort
    autosnip2_AM(opt_autosnip2)
    autosort(snipfile2sortheader(strcat(opt_autosnip2.files{:}(1:end-6),'.ssnp')), datadir);
end
load([datadir '/overview.mat'])
snipdata = sorthead_from_raw(sorthead, opt_autosnip2.snip_channels);
snipdata.unitnames = strrep(cellfun(@(x)['ch',num2str(x)],num2cell(rec_chans),'UniformOutput',false),'.','_');


if manual_clustering
    cass_sort_apply(datadir); %%% this will probably throw error if datadir = working directory
    tempsnips = {}; tempsniptimes = {}; tempunitnames = {};
    for chan = 1:length(snipdata.unitnames)
        if exist([datadir filesep 'chan' num2str(rec_chans(chan))],'dir')  
            tempdir = dir([datadir filesep 'chan' num2str(rec_chans(chan))]);
            [junk junk ext] = cellfun(@fileparts,extractfield(tempdir,'name'),'UniformOutput',false);
            issortfile = strcmp(ext,'.sorted');
            if length(find(issortfile)) > 1
                error('There must be no more than 1 ''.sorted'' file in each chan diretory.')
            elseif length(find(issortfile)) == 1; % if there is one cass sorting file for this chan
                sortedfile = [datadir filesep 'chan' num2str(rec_chans(chan)) filesep tempdir(issortfile).name];
                thisChClusts = sorted_from_raw(merec_obj, sortedfile, opt_autosnip2.snip_range); 
                tempsniptimes = [(flipud(importdata(sortedfile,'-mat')))' tempsniptimes];
                for clust = 1:length(thisChClusts)
                    tempsnips = [thisChClusts(clust) tempsnips]; % preallocate for speed?
%                     tempsniptimes = [importdata(sortedfile,'-mat') tempsniptimes];
                    tempunitnames = [['ch' num2str(rec_chans(chan)) '_cl' num2str(clust)] tempunitnames];
                end
            end
        end
    end
    % Add manually clustered snip data to unclustered snip data. 
    snipdata.snips = [fliplr(tempsnips) snipdata.snips]; 
    snipdata.sniptimes = [fliplr(tempsniptimes) snipdata.sniptimes];
    snipdata.unitnames = [fliplr(tempunitnames) snipdata.unitnames];
    snipdata = rmfield(snipdata,'channels'); % remove to avoid confusion after we add more units
end

par_sets.spikes = NaN(size(par_sets,1),length(snipdata.unitnames));
for event = 1:size(par_sets,1)
    for chan = 1:length(snipdata.unitnames)
        par_sets.spikes(event,chan) = length(find(snipdata.sniptimes{chan} > par_sets.onset(event)...
            & snipdata.sniptimes{chan} <= par_sets.onset(event)+window_scans));
    end
end

%%% Split responses into sf_fixed and tf_fixed, then average the identical
%%% stimulus sets. 
trials_tffixed = par_sets(par_sets.tfreq==pars.tf_fixed,:);
trials_sffixed = par_sets(par_sets.sfreq==pars.sf_fixed,:);
unique_tffixed = size(trials_tffixed,1)/pars.repetitions;
unique_sffixed = size(trials_sffixed,1)/pars.repetitions;
% 'sps_sitemean' contains the spikes for each channel averaged over all
% iterations
% 'sps_shankmean' contains spikes averaged over all sites and iterations
% the 'angleAvg' datasets contain data from the 'tuning' datasets with data
% averaged over all angles
repAvg_tffixed = dataset({NaN(unique_tffixed,5),'Angle','diam','sfreq','tfreq','sps_shankmean'},...
    {cell(unique_tffixed,1),'sps_all'},...
    {NaN(unique_tffixed,length(snipdata.unitnames)), 'sps_sitemean'},...
    'ObsNames',cellstr(num2str((1:unique_tffixed)')));
angleAvg_tffixed = repAvg_tffixed(1:unique_tffixed/pars.nAngles,:); angleAvg_tffixed.Angle=[];
% repAvg_sffixed = dataset({NaN(unique_sffixed,5),'Angle','diam','sfreq','tfreq','sps_shankmean'},...
%     {cell(unique_sffixed,1),'sps_all'},...
%     {NaN(unique_sffixed,length(snipdata.unitnames)), 'sps_sitemean'},...
%     'ObsNames',cellstr(num2str((1:unique_sffixed)')));
% angleAvg_sffixed = repAvg_sffixed(1:unique_sffixed/pars.nAngles,:); angleAvg_sffixed.Angle=[];

comboA = 0; comboB = 0; 
for Angle = 1:pars.nAngles  
    for diam = 1:pars.n_diams
        for sf = 1:pars.n_sfs
            comboA = comboA + 1;
            match=trials_tffixed.sfreq==sf_vals(sf) &...
                  trials_tffixed.diam==diam_vals(diam) &...
                  trials_tffixed.Angle==angle_vals(Angle);
            repAvg_tffixed.sfreq(comboA) = sf_vals(sf);
            repAvg_tffixed.tfreq(comboA) = pars.tf_fixed;
            repAvg_tffixed.diam(comboA) = diam_vals(diam);
            repAvg_tffixed.Angle(comboA) = angle_vals(Angle);
            repAvg_tffixed.sps_all(comboA) = {trials_tffixed.spikes(match,:)};
            repAvg_tffixed.sps_sitemean(comboA,:) = mean(trials_tffixed.spikes(match,:));
            repAvg_tffixed.sps_shankmean(comboA) = mean(repAvg_tffixed.sps_sitemean(comboA,:));
        end
%         for tf = 1:pars.n_tfs
%             comboB = comboB + 1;
%             match=trials_sffixed.tfreq==tf_vals(tf) &...
%                   trials_sffixed.diam==diam_vals(diam) &...
%                   trials_sffixed.Angle==angle_vals(Angle);
%             repAvg_sffixed.sfreq(comboB) = pars.sf_fixed;
%             repAvg_sffixed.tfreq(comboB) = tf_vals(tf);
%             repAvg_sffixed.diam(comboB) = diam_vals(diam);
%             repAvg_sffixed.Angle(comboB) = angle_vals(Angle);
%             repAvg_sffixed.sps_all(comboB) = {trials_sffixed.spikes(match,:)};
%             repAvg_sffixed.sps_sitemean(comboB,:) = mean(trials_sffixed.spikes(match,:));
%             repAvg_sffixed.sps_shankmean(comboB) = mean(repAvg_sffixed.sps_sitemean(comboB,:));
%         end
    end
end

% Get unit-averaged responses for each unit.
anglenames = strrep(cellfun(@(x)['Angle_',num2str(x)],num2cell(angle_vals),'UniformOutput',false),'.','_');
diamnames = strrep(cellfun(@(x)['Diam_',num2str(x)],num2cell(diam_vals),'UniformOutput',false),'.','_'); 
sf_names = strrep(cellfun(@(x)['SF_',num2str(x)],num2cell(sf_vals),'UniformOutput',false),'.','_');
tf_names = strrep(cellfun(@(x)['TF_',num2str(x)],num2cell(tf_vals),'UniformOutput',false),'.','_');
unitsByAngle = repmat({cell(length(snipdata.unitnames),1)},1,length(anglenames)+1);
allAngles_tffixed = dataset(unitsByAngle{:}, NaN(length(snipdata.unitnames),1),...
    'VarNames',[anglenames {'Preferred_Angle_Sps' 'Preferred_Angle'}],'ObsNames', snipdata.unitnames);
% allAngles_sffixed = allAngles_tffixed;

for chan = 1:length(snipdata.unitnames)
    for Angle = 1:pars.nAngles
        allAngles_tffixed{snipdata.unitnames(chan),anglenames(Angle)} =...
            dataset({NaN(pars.n_sfs,pars.n_diams),diamnames{:}},'ObsNames',sf_names);  
%         allAngles_sffixed{snipdata.unitnames(chan),anglenames(Angle)} =...
%             dataset({NaN(pars.n_tfs,pars.n_diams),diamnames{:}},'ObsNames',tf_names);         
        for diam = 1:pars.n_diams
            for sf = 1:pars.n_sfs
                match = repAvg_tffixed.diam == diam_vals(diam) &...
                    repAvg_tffixed.sfreq == sf_vals(sf) &...
                    repAvg_tffixed.Angle == angle_vals(Angle);
                allAngles_tffixed{snipdata.unitnames(chan),anglenames(Angle)}{sf_names(sf),diamnames(diam)} =...
                    repAvg_tffixed.sps_sitemean(match,chan);
            end
%             for tf = 1:pars.n_tfs
%                 match = repAvg_sffixed.diam == diam_vals(diam) &...
%                     repAvg_sffixed.tfreq == tf_vals(tf) &...
%                     repAvg_sffixed.Angle == angle_vals(Angle);
%                 allAngles_sffixed{snipdata.unitnames(chan),anglenames(Angle)}{tf_names(tf),diamnames(diam)} =...
%                     repAvg_sffixed.sps_sitemean(match,chan);
%             end
        end
    end
end

%%%%% find best angle, fill in last 2 rows

% Take average values over all angles. 
comboA = 0; comboB = 0;
for diam = 1:pars.n_diams
    for sf = 1:pars.n_sfs
        comboA = comboA+1;
        match = repAvg_tffixed.diam == diam_vals(diam) &...
                repAvg_tffixed.sfreq == sf_vals(sf); % could limit this match search for speed
        angleAvg_tffixed.diam(comboA) = diam_vals(diam);
        angleAvg_tffixed.sfreq(comboA) = sf_vals(sf);
        angleAvg_tffixed.tfreq(comboA) = pars.tf_fixed;
        angleAvg_tffixed.sps_all(comboA) =... % average the sps_all values for all angles
            {mean(reshape(cell2mat((repAvg_tffixed.sps_all(match))'),...
            pars.repetitions,length(snipdata.unitnames),pars.nAngles),3)};
        angleAvg_tffixed.sps_sitemean(comboA,:) = mean(angleAvg_tffixed.sps_all{comboA,:});
        angleAvg_tffixed.sps_shankmean(comboA) = mean(angleAvg_tffixed.sps_sitemean(comboA,:));
    end
%     for tf = 1:pars.n_tfs
%         comboB = comboB+1;
%         match = repAvg_sffixed.diam == diam_vals(diam) &...
%                 repAvg_sffixed.tfreq == tf_vals(tf); % could limit this match search for speed
%         angleAvg_sffixed.diam(comboB) = diam_vals(diam);
%         angleAvg_sffixed.sfreq(comboB) = pars.sf_fixed;
%         angleAvg_sffixed.tfreq(comboB) = tf_vals(tf);
%         angleAvg_sffixed.sps_all(comboB) =... % average the sps_all values for all angles
%             {mean(reshape(cell2mat((repAvg_sffixed.sps_all(match))'),...
%             pars.repetitions,length(snipdata.unitnames),pars.nAngles),3)};
%         angleAvg_sffixed.sps_sitemean(comboB,:) = mean(angleAvg_sffixed.sps_all{comboB,:});
%         angleAvg_sffixed.sps_shankmean(comboB) = mean(angleAvg_sffixed.sps_sitemean(comboB,:));
%     end
end

% Construct site-specific size tuning curves averaged over angles and repetitions.
unitAvg_tffixed = cell(length(snipdata.unitnames),2);
% unitAvg_sffixed = cell(length(snipdata.unitnames),2);

for chan = 1:length(snipdata.unitnames)
    unitAvg_tffixed{chan} = dataset({NaN(pars.n_sfs,pars.n_diams),diamnames{:}},'ObsNames',sf_names);
%     unitAvg+tffixed{chan,2} = snipdata.unitnames{chan}; 
    unitAvg_sffixed{chan} = dataset({NaN(pars.n_tfs,pars.n_diams),diamnames{:}},'ObsNames',tf_names); 
    for diam = 1:pars.n_diams % make sure sfreqs and tfreqs are taken in the correct order
        for sf = 1:length(sf_vals)
            match = angleAvg_tffixed.diam==diam_vals(diam) & angleAvg_tffixed.sfreq==sf_vals(sf);
            unitAvg_tffixed{chan}.(diamnames{diam})(sf_names{sf}) = angleAvg_tffixed.sps_sitemean(match,chan);
        end
%         for tf = 1:length(tf_vals)
%             match = angleAvg_sffixed.diam==diam_vals(diam) & angleAvg_sffixed.tfreq==tf_vals(tf);
%             unitAvg_sffixed{chan}.(diamnames{diam})(tf_names{tf}) = angleAvg_sffixed.sps_sitemean(match,chan);
%         end
    end
end
unitAvg_tffixed(:,2) = snipdata.unitnames'; % add channel labels 
% unitAvg_sffixed(:,2) = snipdata.unitnames'; % add channel labels 


% Construct size tuning curves averaged over angles, repetitions, and shanks. 
grandAvg_tffixed = dataset({NaN(pars.n_sfs,pars.n_diams),diamnames{:}},'ObsNames',sf_names);
% grandAvg_sffixed = dataset({NaN(pars.n_tfs,pars.n_diams),diamnames{:}},'ObsNames',tf_names); 

for diam = 1:pars.n_diams % make sure sfreqs and tfreqs are taken in the correct order
    for sf = 1:length(sf_vals)
        match = angleAvg_tffixed.diam==diam_vals(diam) & angleAvg_tffixed.sfreq==sf_vals(sf);
        grandAvg_tffixed.(diamnames{diam})(sf_names{sf}) = angleAvg_tffixed.sps_shankmean(match);
    end
%     for tf = 1:length(tf_vals)
%         match = angleAvg_sffixed.diam==diam_vals(diam) & angleAvg_sffixed.tfreq==tf_vals(tf);
%         grandAvg_sffixed.(diamnames{diam})(tf_names{tf}) = angleAvg_sffixed.sps_shankmean(match);
%     end
end

[junk stim_filename junk] = fileparts(stim_file);
save([stim_filename '_analyzed.mat']);

%% Plotting
% Plot size tuning curves for fixed tf for one unit
% % % % % plotsite = '13'; % number of site in merec indices as string
% % % % % for i = 1:length(anglenames)
% % % % %     figure; plot(transpose(double(allAngles_tffixed{['Ch_' plotsite],anglenames(i)}))/pars.stimdur);
% % % % %     set(gca,'XTickLabel',diamnames);
% % % % %     legend(sf_names,'Interpreter','none')
% % % % %     title(['ch' plotsite ' ' anglenames(i)],'Interpreter','none')
% % % % %     ylabel('hz')
% % % % % end

end