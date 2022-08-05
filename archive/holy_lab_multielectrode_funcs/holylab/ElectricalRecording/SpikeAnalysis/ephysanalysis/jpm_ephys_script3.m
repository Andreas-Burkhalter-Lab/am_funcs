% jpm_ephys_script3.m: 2008_01_24 script to extract and plot a running tab
% of delta_r values from the AOB sulfated steroid project (2008_01_24 -
% ***)

%% Load data into ephys (from generic epanalyze1.m)
% clear any currently loaded variables from memory/figures
clear;
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
%    data.valvelabels(1) = {'50 mM KCl'}; data.valvelabels(16) = []; % this expt had an error in the user comment section
%    [data.valvelabels] = deal(jpm_simple_mixed_urine_valvelabels_2007_04_01);  % UNCOMMENT FOR MANUAL VALVE LABELING
    data.valvelabels(1) = {'Ringer''s'};  % mistake on 2008_03_06
    for valve_idx = 1:size(data(1).valvelabels,2)
        if isempty(data(1).valvelabels{valve_idx})
            for data_idx = 1:size(data,2)
                data(data_idx).valvelabels(valve_idx) = [];
            end
        end
    end
    valvelabels = {data(1).valvelabels};
     
% set the names for the subdirectories containing sorted data information
    [data.sortfile] = deal('sort1');
    chanfile.channel = 1; % the channel number to sort
    chanfile.file = '1'; % the basname of the .sorted file to load

% update data using 'ephys_cass_fetch'
    data = ephys_cass_fetch(data, chanfile);
    [data.sort_cass_chanfile] = deal(chanfile);
    data = ephysfetch(data, 'celltimes');
    
% 'trange' contains a 2-vector (in seconds) which delimits the time 
% before and after your stimulus onset to consider for plotting:
    trange = [-10 20];

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
    
%%  NOW, use this to extract delta_r and selectivity the site (sniptimes or celltimes only)

% extract delta_r information on each valve for this experiment
% make decision to use 'celltimes' or 'sniptimes'
choice = 'celltimes';
if choice == 'celltimes'
    number_of_cells = size(data(1).cellnums,2);
    cellnums = data(1).cellnums;
    for idx_label = 1:size(valvelabels{1},2)
        temp_delta_r = calc_delta_r_from_ephys(data, valvelabels{1}{idx_label}, [1 7], 'celltimes');
        if ~isempty(temp_delta_r)
            for cell_idx = 1:number_of_cells
                delta_r{cell_idx}{idx_label} = temp_delta_r{cell_idx};
                delta_r_mean{cell_idx}(idx_label) = mean(delta_r{cell_idx}{idx_label});
                delta_r_stderr{cell_idx}(idx_label) = (std(delta_r{cell_idx}{idx_label}))/(sqrt(length(delta_r{cell_idx}{idx_label})));
            end
        else
            for cell_idx = 1:number_of_cells
                delta_r{cell_idx}{idx_label} = NaN;
                delta_r_mean{cell_idx}(idx_label) = NaN;
                delta_r_stderr{cell_idx}(idx_label) = NaN;
            end
        end
    end
elseif choice == 'sniptimes'
    number_of_cells = NaN;
    cellnums = NaN;
    for idx_label = 1:size(valvelabels{1},2)
        delta_r(idx_label) = {calc_delta_r_from_ephys(data, valvelabels{1}{idx_label}, [1 7], 'sniptimes')};
        delta_r_mean(idx_label) = mean(delta_r{idx_label});
        delta_r_stderr(idx_label) = (std(delta_r{idx_label}))/(sqrt(length(delta_r{idx_label})));
    end
end

% calculate selectivity index (female-positive) @ 1:100 concentration
% grab the indices we want (a little laborious considering the job is done in utags, identities)
for idx = 1:size(valvelabels{1},2)
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
if ~isnan(number_of_cells)
    for cell_idx = 1:number_of_cells
        fem_selectivity(cell_idx) = (delta_r_mean{cell_idx}(fem_index)-delta_r_mean{cell_idx}(mal_index))...
                  /(abs(delta_r_mean{cell_idx}(fem_index))+abs(delta_r_mean{cell_idx}(mal_index)));
    end
else          
    fem_selectivity = (delta_r_mean(fem_index)-delta_r_mean(mal_index))...
                 /(abs(delta_r_mean(fem_index))+abs(delta_r_mean(mal_index)));
end


%% calculate lifetime sparseness (in progress)
% define single compounds used
% NOT USED AT MOMENT
  %sulf_steroid_list = {'A6940','A7010','A7864',...
  %                   'E0893','E1050','E4105',...
  %                   'P3317','P3865','P8200',...
  %                   'Q1570','Q3383','Q3910'};

%identify valves with sulfated steroids (using Steraloids code strategy)              
sulf_steroids_used = [];
for idx = 1:size(valvelabels{1},2)
    if ~isempty(valvelabels{1}{idx})
        if~isempty(regexp(valvelabels{1}{idx},'[A E P Q]\d\d\d\d'))
            if isempty(sulf_steroids_used)
                sulf_steroids_used = idx;
            else
                sulf_steroids_used(end+1) = idx;
            end
        end
    end
end
% 
% % if multiple cells, this means we need to parse the data to make it more
% % easily useable (change delta_r from {valve}->(cell) to {cell}->(valves)
% for valve_idx = 1:size(valvelabels{1},2)
%     for cell_idx = 1:number_of_cells
%         conv_delta_r_mean{cell_idx}(valve_idx) = delta_r_mean{cell_idx}{valve_idx};
%     end
% end

if ~isnan(number_of_cells)
    for cell_idx = 1:number_of_cells
        sulf_steroid_lifetime_sparseness(cell_idx) = ...
         (1 - (sum(delta_r_mean{cell_idx}(sulf_steroids_used)/length(sulf_steroids_used))^2 ...
                /sum((delta_r_mean{cell_idx}(sulf_steroids_used).^2)/length(sulf_steroids_used)))) ...
                /(1-(1/length(sulf_steroids_used)));
    end
else
    sulf_steroid_lifetime_sparseness = ...
        (1 - (sum(delta_r_mean{1}(sulf_steroids_used)/length(sulf_steroids_used))^2 ...
        /sum((delta_r_mean{1}(sulf_steroids_used).^2)/length(sulf_steroids_used)))) ...
        /(1-(1/length(sulf_steroids_used)));
end
   
%% calculate statistics (Student's 2-tailed, unpaired t-test; ranksum test)
for idx = 1:size(valvelabels{1},2)
    if exist('ring_index','var')
        if ~isnan(number_of_cells)
            for cell_idx = 1:number_of_cells
                if ~isempty(delta_r{cell_idx}{idx})&&length(delta_r{cell_idx}{idx})>=3 ...
                        &&length(delta_r{cell_idx}{ring_index})>=3
                    [h, t_test_p{cell_idx}(idx)] = ttest2(delta_r{cell_idx}{idx}, delta_r{cell_idx}{ring_index}, .05, 'both');
                    ranksum_p{cell_idx}(idx) = ranksum(delta_r{cell_idx}{idx},delta_r{cell_idx}{ring_index});
                else
                    t_test_p{cell_idx}(idx) = NaN;
                    ranksum_p{cell_idx}(idx) = NaN;
                end
            end
        else   %% PROBABLY BROKEN RIGHT NOW!!!!
            if ~isempty(delta_r{idx})&&length(delta_r{idx})>=3 ...
                        &&length(delta_r{ring_index})>=3
                    [h, t_test_p{cell_idx}(idx)] = ttest2(delta_r{idx}(cell_idx), delta_r{ring_index}(cell_idx), .05, 'both');
                    ranksum_p{cell_idx}(idx) = ranksum(delta_r{idx}(cell_idx),delta_r{ring_index}(cell_idx));
                else
                    t_test_p{cell_idx}(idx) = NaN;
                    ranksum_p{cell_idx}(idx) = NaN;
            end
        end
    end
end

if ~isnan(number_of_cells)
    for cell_idx = 1:number_of_cells
        if sum(t_test_p{cell_idx}) == 0; t_test_p{cell_idx} = NaN;end
        if sum(ranksum_p{cell_idx}) == 0; ranksum_p{cell_idx} = NaN;end
    end
else
    
    if sum(t_test_p) == 0; t_test_p = {NaN};end
    if sum(ranksum_p) == 0; ranksum_p = {NaN};end
end
    
%% put together save structure
analysis = struct;
if ~isnan(number_of_cells)
    for cell_idx = 1:number_of_cells
        analysis(cell_idx).fem_selectivity_dilution = fem_selectivity_dilution;
        analysis(cell_idx).fem_selectivity = fem_selectivity(cell_idx);
        analysis(cell_idx).delta_r = delta_r{cell_idx};
        analysis(cell_idx).valvelabels = valvelabels{1};
        analysis(cell_idx).delta_r_mean = delta_r_mean{cell_idx};
        analysis(cell_idx).delta_r_stderr = delta_r_stderr{cell_idx};
        analysis(cell_idx).sulf_steroid_lifetime_sparseness = sulf_steroid_lifetime_sparseness(cell_idx);
        analysis(cell_idx).sulf_steroid_valves = sulf_steroids_used;
        analysis(cell_idx).t_test_pvalue = t_test_p{cell_idx};
        analysis(cell_idx).ranksum_p = ranksum_p{cell_idx};
        analysis(cell_idx).cellnum = cellnums(cell_idx);
    end
else
    analysis.fem_selectivity_dilution = fem_selectivity_dilution;
    analysis.fem_selectivity = fem_selectivity;
    analysis.delta_r = delta_r;
    analysis.valvelabels = valvelabels{1};
    analysis.delta_r_mean = delta_r_mean;
    analysis.delta_r_stderr = delta_r_stderr;
    analysis.sulf_steroid_lifetime_sparseness = sulf_steroid_lifetime_sparseness;
    analysis.sulf_steroid_valves = sulf_steroids_used;
    analysis.t_test_pvalue = t_test_p;
    analysis.ranksum_p = ranksum_p;
end

%% save data into directory
save('analysis.mat', 'analysis');