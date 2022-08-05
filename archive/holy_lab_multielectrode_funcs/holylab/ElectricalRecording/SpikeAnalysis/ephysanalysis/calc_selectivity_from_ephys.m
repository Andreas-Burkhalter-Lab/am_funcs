% calc_selectivity_from_ephys (JPM, script only at this time)

%% Load data into ephys (from generic epanalyze1.m)
% clear any currently loaded variables from memory/figures
clear;
clf;
% gather names of the cleaned snippet files in the current directory
    files = dirbyname('*.ssnp');
% gather names of the .vlv files containing stimulus timing info from the
% current directory
    stimfilenames = dirbyname('*.vlv');
% in case naming convention does not load files in time order, force the
% issue:
    files = sort_ai_by_time(files);
% initialize the ephys structure for the files in this directory using the
% file header:
    data = ephysfromai(files, struct('usefullpath',false));
% update the current 'data' (i.e. the current ephys struct) with the
% stimulus information (timing, valve numbers)
    data = ephysfetch(data,'stimulus');
    data = ephysfetch(data,'sniptimes');
    
% create a cell array of the available valve labels from the ephys
% structure:
    data.valvelabels(1) = {'50 mM KCl'}; 
    valvelabels = {data.valvelabels}; % this expt had an error in the user comment section
     
% set the names for the subdirectories containing sorted data information
    [data.sortfile] = deal('sort1');
    chanfile.channel = 1; % the channel number to sort
    chanfile.file = '1'; % the basname of the .sorted file to load

% update data using 'ephys_cass_fetch'
    data = ephys_cass_fetch(data, chanfile);
    [data.sort_cass_chanfile] = deal(chanfile);
    data = ephysfetch(data, 'celltimes');
    %data = ephysfetch(data, 'cellnums');
    
% 'trange' contains a 2-vector (in seconds) which delimits the time 
% before and after your stimulus onset to consider for plotting:
    trange = [-10 25];

% call 'intervalsfromstim' to assign arrays containing the intervals,
% valve identities, and valve names from this experiment
    [intervals,identities,vlvsout] = intervalsfromstim(data,trange);

% call ephyssubrange to update your data (ephys struct) to only
% contain the sniptimes, etc from your designated intervals
    data = ephyssubrange(data,intervals);

% one-line 'if' statement to collapse your data across cells into a
% one-dimensional cell array:
    if iscell(data), data = cat(2,data{:}); end

% call 'deal' (MATLAB function) to add the offset time from trange
% to every cell in the data array (this allows the time of stimulus
% onset to be equal to zero)
    [data.toffset] = deal(trange(1));

% collapes the valve identities corresponding to each
% interval into a 1-dimensional cell array
    identities = cat(2,identities{:});

% uncomment the following lines to use stimulus timing in the tags
    %tagops.addtime = 1;         % flag to add the time
    %tagops.timedelim = ' $ ';   % will occur on plot before numeral
    %tagops.timeprec = 0.5;      % 0.5 sec precision (sloppy)
    
% set the "tag" of each interval to the name of the valve label
% 'utags' will become a cell array of strings from all used valves
    %[data,utags] = ephystag(data,valvelabels{1}(identities), tagops);
    [data,utags] = ephystag(data,valvelabels{1}(identities));
    
%% Now calculate the selectivity index (female-based)

% choose 1:100 static for the moment

valvelabels{1}(16)=[];  % for some reason this has an empty char @ end (may be bug)

% MAY BE A USEFUL FUNCTION SOON: grab indices based on valvelabel name
% (esp. once names have been normalized using the .xdb)
fem100intervals = nonzeros(find(strmatch('1:100 Balb FU',valvelabels{1}(:),'exact')==identities));
mal100intervals = nonzeros(find(strmatch('1:100 Balb MU',valvelabels{1}(:),'exact')==identities));

% extract the data for just those intervals you are interested in:
fem100data = data(fem100intervals);
mal100data = data(mal100intervals);

% choose the range (in seconds) you wish to integrate) post-stimulus
% example: [1 12] gets 11 seconds worth of spikes starting 1 second after the stimulus
delta_r_range = [1 11];

% calculate the delta r value for fem and mal100
for idx_f = 1:size(fem100data,2);
    delta_r_scanrange_f{idx_f} = (delta_r_range-fem100data(idx_f).toffset)*fem100data(idx_f).scanrate+fem100data(idx_f).scanrange(1);
    base_scanrange_f{idx_f} = [fem100data(idx_f).scanrange(1) 
                               fem100data(idx_f).scanrange(1)+...
                                     (delta_r_range(2)-delta_r_range(1))*fem100data(idx_f).scanrate];
    numfem100spikes(idx_f) = length(between(fem100data(idx_f).sniptimes{:},...
                             delta_r_scanrange_f{idx_f}(1),...
                             delta_r_scanrange_f{idx_f}(2)));
    basefem100spikes(idx_f) = length(between(fem100data(idx_f).sniptimes{:},...
                             base_scanrange_f{idx_f}(1),...
                             base_scanrange_f{idx_f}(2)));
    stim_firerate_f = numfem100spikes(idx_f)/(delta_r_range(2)-delta_r_range(1));
    base_firerate_f = basefem100spikes(idx_f)/(delta_r_range(2)-delta_r_range(1));
    delta_r_f(idx_f) = stim_firerate_f-base_firerate_f;
end

% calculate the delta r value for fem and mal100
for idx_m = 1:size(mal100data,2);
    delta_r_scanrange_m{idx_m} = (delta_r_range-mal100data(idx_m).toffset)*mal100data(idx_m).scanrate+mal100data(idx_m).scanrange(1);
    base_scanrange_m{idx_m} = [mal100data(idx_m).scanrange(1) 
                               mal100data(idx_m).scanrange(1)+...
                                     (delta_r_range(2)-delta_r_range(1))*mal100data(idx_m).scanrate];
    nummal100spikes(idx_m) = length(between(mal100data(idx_m).sniptimes{:},...
                             delta_r_scanrange_m{idx_m}(1),...
                             delta_r_scanrange_m{idx_m}(2)));
    basemal100spikes(idx_m) = length(between(mal100data(idx_m).sniptimes{:},...
                             base_scanrange_m{idx_m}(1),...
                             base_scanrange_m{idx_m}(2)));
    stim_firerate_m = nummal100spikes(idx_m)/(delta_r_range(2)-delta_r_range(1));
    base_firerate_m = basemal100spikes(idx_m)/(delta_r_range(2)-delta_r_range(1));
    delta_r_m(idx_m) = stim_firerate_m-base_firerate_m;
end

% set up a struct for saving
save_data = struct;
save_data.delta_r = struct;
save_data.basefilename = data(1).basefilename;
save_data.analyzedate = datenum(now);
save_data.delta_r.MU100 = delta_r_m;
save_data.delta_r.FU100 = delta_r_f;

% now that the data is saved, calculate avg and SEM for this cell
save_data.delta_r_avg.MU100 = mean(delta_r_m);
save_data.delta_r_avg.FU100 = mean(delta_r_f);
save_data.delta_r_sem.MU100 = std(delta_r_m)/sqrt(length(delta_r_m));
save_data.delta_r_sem.FU100 = std(delta_r_f)/sqrt(length(delta_r_f));

save_data.selectivity = (save_data.delta_r_avg.FU100-save_data.delta_r_avg.MU100)/...
                         (abs(save_data.delta_r_avg.FU100)+abs(save_data.delta_r_avg.MU100));

save('analysis.mat','save_data');
