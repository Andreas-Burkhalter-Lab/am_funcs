% jpm_ephys_script2.m: 11-30-07 script to extract and plot the things I
% want for the NRSA Dec 2007 deadline

%% Load data into ephys (from generic epanalyze1.m)
% clear any currently loaded variables from memory/figures
clear;
clf;
% gather names of the cleaned snippet files in the current directory
    %files = dirbyname('*clean.ssnp');
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
%    data.valvelabels(1) = {'50 mM KCl'}; % this expt had an error in the user comment section
    [data.valvelabels] = deal(jpm_simple_mixed_urine_valvelabels_2007_03_06);
    valvelabels = {data.valvelabels};
     
% set the names for the subdirectories containing sorted data information
    [data.sortfile] = deal('sort1');
    chanfile.channel = 1; % the channel number to sort
    chanfile.file = '1'; % the basname of the .sorted file to load

% update data using 'ephys_cass_fetch'
    %data = ephys_cass_fetch(data, chanfile);
    %[data.sort_cass_chanfile] = deal(chanfile);
    %data = ephysfetch(data, 'celltimes');
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
    
%%  NOW, use this to extract delta_r and selectivity the site (snippets only)

% extract delta_r information on each valve for this experiment
for idx_label = 1:size(valvelabels{1},2)
    delta_r(idx_label) = {calc_delta_r_from_ephys(data, valvelabels{1}{idx_label}, [1 11])};
    delta_r_mean(idx_label) = mean(delta_r{idx_label});
    delta_r_stderr(idx_label) = (std(delta_r{idx_label}))/(sqrt(length(delta_r{idx_label})));
end

% calculate selectivity index (female-positive) @ 1:100 concentration
% grab the indices we want (a little laborious considering the job is done in utags, identities)
for idx = 1:size(valvelabels{1},2);
    if ~isempty(valvelabels{1}{idx})
        fem = strmatch('1:100 Balb FU', valvelabels{1}{idx});
        mal = strmatch('1:100 Balb MU', valvelabels{1}{idx});
        ring = strmatch('Ringer''s', valvelabels{1}{idx});
        if ~isempty(fem)
            if ~exist('fem_index','var')
                fem_index = idx;
            end
        end
        if ~isempty(mal)
            if ~exist('mal_index', 'var')
                mal_index = idx;
            end
        end
        if ~isempty(ring)
            if ~exist('ring_index', 'var')
                ring_index = idx;
            end
        end
    end
end

%% calculate selectivity
fem_selectivity_dilution = 100;
fem_selectivity = (delta_r_mean(fem_index)-delta_r_mean(mal_index))...
                  /(abs(delta_r_mean(fem_index))+abs(delta_r_mean(mal_index)));

%% calculate statistics (Student's 2-tailed, unpaired t-test)
t_test_p = zeros(size(valvelabels{1}));

for idx = 1:size(valvelabels{1},2)
    if exist('ring_index','var')
        if ~isempty(delta_r{idx})&&length(delta_r{idx})==length(delta_r{ring_index})
            [h, t_test_p(idx)] = ttest2(delta_r{idx}, delta_r{ring_index}, .05, 'both');
        else
            t_test_p(idx) = NaN;
        end
    elseif ~isempty(delta_r{idx})
        [h, t_test_p(idx)] = ttest(delta_r{idx});
    end
end
%% put together save structure
analysis = struct;
analysis.fem_selectivity_dilution = fem_selectivity_dilution;
analysis.fem_selectivity = fem_selectivity;
analysis.delta_r = delta_r;
analysis.valvelabels = valvelabels{1};
analysis.delta_r_mean = delta_r_mean;
analysis.delta_r_stderr = delta_r_stderr;
analysis.t_test_pvalue = t_test_p;

%% save data into directory
save('analysis.mat', 'analysis');