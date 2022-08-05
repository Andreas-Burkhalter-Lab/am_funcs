function fmb_analysis = analyze_fmb(directory, conditions, condition_code)
% ANALYZE_FMB takes .fmbn files in a directory and creates a 'fmb_analysis' structure
%  *** very much a work in progress***
% Syntax: analyze_fmb(directory <string>,...
%                     conditions <cell array of 2 strings>
%                     condition_code <2 x 2 matrix>);
%     'directory' is the directory containing your .fmbn files
%     'conditions' is a descriptive (2-element) cell array of strings
%         corresponding to the 'a' and 'b' conditions in the fmb data
%     'condition_code' is a 4x2x2 matrix of integers containing the
%         controller number.  
%         use:(tricky): 
%            condition_code(1,:,1) contains the [x y] integer labels of axes in condition 'a' in file 1, 
%            condition_code(2,:,1) contains ''                               in condition 'a' in file 2,
%            condition_code(3,:,2) contains ''                               in condition 'b' in file 3,
%            and so on...

% Copyright 2008 Julian Meeks (Timothy Holy Laboratory)       
%
% Revision History:
% 2008_06_17: Wrote it (JPM)

if exist([directory '/fmb_analysis.mat'], 'file')
    load([directory '/fmb_analysis.mat'], '-mat');
else
    fmb_analysis = struct;
end

% Version #
fmba_version = 1.0;

if iscellstr(conditions)
    fmb_analysis.conditions = {'1:2 Balb FU', 'Ringer''s nosugar'};
else
    fmb_analysis.conditions = {'', ''};
end

filenames = dirbyname([directory '/*.fmbn']);

% currently only supports when all experiments have the same 'key'  in
% future will need to be the same size as 'filenames' below (if 'key' to
% conditions varies across experiments
if size(condition_code, 1) ~= size(filenames,2)
    error('''condition_code'' mismatch.');
end


for f_idx = 1:size(filenames,2);
    load(filenames{f_idx}, '-mat');
    fmb_analysis.filenames{f_idx} = filenames{f_idx};
    
    % determine raw waveforms divided by minimum difference value (337)
    fmb_analysis.raw{f_idx} = flipdim(rot90(fmb.combined_nonsparse_events./337),1);
    
    % determine per-controller difference and sum(abs(diff())) values ('diff' and 'int_diff')
    fmb_analysis.diff{f_idx} = diff(fmb_analysis.raw{f_idx},1);
    count = 1;
    for i_idx = 101:100:floor(size(fmb_analysis.diff{f_idx},1))-101
        fmb_analysis.int_diff_1s{f_idx}(count,:) = ...
            sum(abs(fmb_analysis.diff{f_idx}(i_idx-100:i_idx+99,:)),1);
        fmb_analysis.int_diff_1s{f_idx}(count,:) = ...
            sum(abs(fmb_analysis.diff{f_idx}(i_idx-100:i_idx+99,:)),1);
        count = count+1;
    end
    count = 1;
    for i_idx = 1001:1000:floor(size(fmb_analysis.diff{f_idx},1))-1001
        fmb_analysis.int_diff_10s{f_idx}(count,:) = ...
            sum(abs(fmb_analysis.diff{f_idx}(i_idx-1000:i_idx+999,:)),1);
        fmb_analysis.int_diff_10s{f_idx}(count,:) = ...
            sum(abs(fmb_analysis.diff{f_idx}(i_idx-1000:i_idx+999,:)),1);
        count = count+1;
    end
    
    % determine the per-condition difference and 'int_diff'
    fmb_analysis.combined_diff{f_idx}(:,1) = ...
        sum(fmb_analysis.diff{f_idx}(:, condition_code(f_idx,:,1)),2);
    fmb_analysis.combined_diff{f_idx}(:,2) = ...
        sum(fmb_analysis.diff{f_idx}(:, condition_code(f_idx,:,2)),2);
    fmb_analysis.combined_int_diff_1s{f_idx}(:,1) =...
        sum(abs(fmb_analysis.int_diff_1s{f_idx}(:, condition_code(f_idx,:,1))),2);
    fmb_analysis.combined_int_diff_1s{f_idx}(:,2) =...
        sum(abs(fmb_analysis.int_diff_1s{f_idx}(:, condition_code(f_idx,:,2))),2);
    fmb_analysis.combined_int_diff_10s{f_idx}(:,1) =...
        sum(abs(fmb_analysis.int_diff_10s{f_idx}(:, condition_code(f_idx,:,1))),2);
    fmb_analysis.combined_int_diff_10s{f_idx}(:,2) =...
        sum(abs(fmb_analysis.int_diff_10s{f_idx}(:, condition_code(f_idx,:,2))),2);
    
    % approximate the number of large events (>3x stdev) and small events (stdev<x<3xstdev)
    % during the recording (somewhat arbitrary)
    
    % per-controller
    for idx = 1:size(fmb_analysis.int_diff_1s{f_idx},2)
        % for 1s integration windows:
        stdev = std(fmb_analysis.int_diff_1s{f_idx}(:,idx));
        fmb_analysis.approx_events.large_1s(f_idx,idx) =...
            size(find(fmb_analysis.int_diff_1s{f_idx}(:,idx) > 3*stdev),1);
        fmb_analysis.approx_events.small_1s(f_idx,idx) =...
            size(...
               intersect(...
                  find(fmb_analysis.int_diff_1s{f_idx}(:,idx) < 3*stdev),...
                  find(fmb_analysis.int_diff_1s{f_idx}(:,idx) > stdev)...
               )...
            ,1);
        
        fmb_analysis.approx_events.total_1s(f_idx, :) = ...
            fmb_analysis.approx_events.large_1s(f_idx,:) + fmb_analysis.approx_events.small_1s(f_idx,:);
        % for 10s integration windows:
        stdev = std(fmb_analysis.int_diff_10s{f_idx}(:,idx));
        fmb_analysis.approx_events.large_10s(f_idx,idx) =...
            size(find(fmb_analysis.int_diff_10s{f_idx}(:,idx) > 3*stdev),1);
        fmb_analysis.approx_events.small_10s(f_idx,idx) =...
            size(...
               intersect(...
                  find(fmb_analysis.int_diff_10s{f_idx}(:,idx) < 3*stdev),...
                  find(fmb_analysis.int_diff_10s{f_idx}(:,idx) > stdev)...
               )...
            ,1);
        % combine large and small events:
        fmb_analysis.approx_events.total_10s(f_idx, :) = ...
            fmb_analysis.approx_events.large_10s(f_idx,:) + fmb_analysis.approx_events.small_10s(f_idx,:);
    end
    
    % per-condition
    for idx = 1:size(fmb_analysis.combined_int_diff_1s{f_idx},2)
        % for 1s integration windows:
        stdev = std(fmb_analysis.combined_int_diff_1s{f_idx}(:,idx));
        fmb_analysis.approx_events.combined_large_1s(f_idx,idx) =...
            size(find(fmb_analysis.combined_int_diff_1s{f_idx}(:,idx) > 3*stdev),1);
        fmb_analysis.approx_events.combined_small_1s(f_idx,idx) =...
            size(...
               intersect(...
                  find(fmb_analysis.combined_int_diff_1s{f_idx}(:,idx) < 3*stdev),...
                  find(fmb_analysis.combined_int_diff_1s{f_idx}(:,idx) > stdev)...
               )...
            ,1);
        % combine large and small events:
        fmb_analysis.approx_events.combined_total_1s(f_idx, :) = ...
            fmb_analysis.approx_events.combined_large_1s(f_idx,:) + fmb_analysis.approx_events.combined_small_1s(f_idx,:);
        % for 10s integration windows:
        stdev = std(fmb_analysis.combined_int_diff_10s{f_idx}(:,idx));
        fmb_analysis.approx_events.combined_large_10s(f_idx,idx) =...
            size(find(fmb_analysis.combined_int_diff_10s{f_idx}(:,idx) > 3*stdev),1);
        fmb_analysis.approx_events.combined_small_10s(f_idx,idx) =...
            size(...
               intersect(...
                  find(fmb_analysis.combined_int_diff_10s{f_idx}(:,idx) < 3*stdev),...
                  find(fmb_analysis.combined_int_diff_10s{f_idx}(:,idx) > stdev)...
               )...
            ,1);
        fmb_analysis.approx_events.combined_total_10s(f_idx, :) = ...
            fmb_analysis.approx_events.combined_large_10s(f_idx,:) + fmb_analysis.approx_events.combined_small_10s(f_idx,:);

    end
    
    fmb_analysis.sum_activity_1s(f_idx,:) = sum(fmb_analysis.int_diff_1s{f_idx}, 1);
    fmb_analysis.sum_activity_10s(f_idx,:) = sum(fmb_analysis.int_diff_10s{f_idx}, 1);
    fmb_analysis.mean_activity_1s(f_idx,:) = mean(fmb_analysis.int_diff_1s{f_idx}, 1);
    fmb_analysis.mean_activity_10s(f_idx,:) = mean(fmb_analysis.int_diff_10s{f_idx}, 1);
    fmb_analysis.peak_activity_1s(f_idx,:) = max(fmb_analysis.int_diff_1s{f_idx});
    fmb_analysis.peak_activity_10s(f_idx,:) = max(fmb_analysis.int_diff_10s{f_idx});

    fmb_analysis.combined_sum_activity_1s(f_idx,:) = sum(fmb_analysis.combined_int_diff_1s{f_idx}, 1);
    fmb_analysis.combined_mean_activity_1s(f_idx,:) = mean(fmb_analysis.combined_int_diff_1s{f_idx}, 1);
    fmb_analysis.combined_mean_activity_10s(f_idx,:) = mean(fmb_analysis.combined_int_diff_10s{f_idx}, 1);
    fmb_analysis.combined_peak_activity_1s(f_idx,:) = max(fmb_analysis.combined_int_diff_1s{f_idx});
    fmb_analysis.combined_peak_activity_10s(f_idx,:) = max(fmb_analysis.combined_int_diff_10s{f_idx});

end


% do stats! paired, two-tailed t-tests for now
a = cat(1,diag(fmb_analysis.sum_activity_1s(:,condition_code(:,1,1))), diag(fmb_analysis.sum_activity_1s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.sum_activity_1s(:,condition_code(:,1,2))), diag(fmb_analysis.sum_activity_1s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.sum_activity_1s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.sum_activity_1s = p;
a = cat(1,diag(fmb_analysis.sum_activity_10s(:,condition_code(:,1,1))), diag(fmb_analysis.sum_activity_10s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.sum_activity_10s(:,condition_code(:,1,2))), diag(fmb_analysis.sum_activity_10s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.sum_activity_10s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.sum_activity_10s = p;

a = cat(1,diag(fmb_analysis.mean_activity_1s(:,condition_code(:,1,1))), diag(fmb_analysis.mean_activity_1s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.mean_activity_1s(:,condition_code(:,1,2))), diag(fmb_analysis.mean_activity_1s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.mean_activity_1s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.mean_activity_1s = p;
a = cat(1,diag(fmb_analysis.mean_activity_10s(:,condition_code(:,1,1))), diag(fmb_analysis.mean_activity_10s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.mean_activity_10s(:,condition_code(:,1,2))), diag(fmb_analysis.mean_activity_10s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.mean_activity_10s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.mean_activity_10s = p;

a = cat(1,diag(fmb_analysis.peak_activity_1s(:,condition_code(:,1,1))), diag(fmb_analysis.peak_activity_1s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.peak_activity_1s(:,condition_code(:,1,2))), diag(fmb_analysis.peak_activity_1s(:,condition_code(:,2,2))));
[h p] = ttest(a, b,0.05, 'both');
fmb_analysis.ttest.peak_activity_1s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.peak_activity_1s = p;
a = cat(1,diag(fmb_analysis.peak_activity_10s(:,condition_code(:,1,1))), diag(fmb_analysis.peak_activity_10s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.peak_activity_10s(:,condition_code(:,1,2))), diag(fmb_analysis.peak_activity_10s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.peak_activity_10s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.peak_activity_10s = p;

a = cat(1,diag(fmb_analysis.approx_events.large_1s(:,condition_code(:,1,1))), diag(fmb_analysis.approx_events.large_1s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.approx_events.large_1s(:,condition_code(:,1,2))), diag(fmb_analysis.approx_events.large_1s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.approx_events_large_1s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.approx_events_large_1s = p;
a = cat(1,diag(fmb_analysis.approx_events.large_10s(:,condition_code(:,1,1))), diag(fmb_analysis.approx_events.large_10s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.approx_events.large_10s(:,condition_code(:,1,2))), diag(fmb_analysis.approx_events.large_10s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.approx_events_large_10s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.approx_events_large_10s = p;

a = cat(1,diag(fmb_analysis.approx_events.small_1s(:,condition_code(:,1,1))), diag(fmb_analysis.approx_events.small_1s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.approx_events.small_1s(:,condition_code(:,1,2))), diag(fmb_analysis.approx_events.small_1s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.approx_events_small_1s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.approx_events_small_1s = p;
a = cat(1,diag(fmb_analysis.approx_events.small_10s(:,condition_code(:,1,1))), diag(fmb_analysis.approx_events.small_10s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.approx_events.small_10s(:,condition_code(:,1,2))), diag(fmb_analysis.approx_events.small_10s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.approx_events_small_10s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.approx_events_small_10s = p;

a = cat(1,diag(fmb_analysis.approx_events.total_1s(:,condition_code(:,1,1))), diag(fmb_analysis.approx_events.total_1s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.approx_events.total_1s(:,condition_code(:,1,2))), diag(fmb_analysis.approx_events.total_1s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.approx_events_total_1s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.approx_events_total_1s = p;
a = cat(1,diag(fmb_analysis.approx_events.total_10s(:,condition_code(:,1,1))), diag(fmb_analysis.approx_events.total_10s(:,condition_code(:,2,1))));
b = cat(1,diag(fmb_analysis.approx_events.total_10s(:,condition_code(:,1,2))), diag(fmb_analysis.approx_events.total_10s(:,condition_code(:,2,2))));
[h p] = ttest(a, b, 0.05, 'both');
fmb_analysis.ttest.approx_events_total_10s = p;
[h p] = ttest2(a, b, 0.05, 'both');
fmb_analysis.ttest2.approx_events_total_10s = p;

% if hand_scored... exist, do the t_tests on them too
if isfield(fmb_analysis, 'hand_scored_events')
    a = cat(1,diag(fmb_analysis.hand_scored_events(:,condition_code(:,1,1))), diag(fmb_analysis.hand_scored_events(:,condition_code(:,2,1))));
    b = cat(1,diag(fmb_analysis.hand_scored_events(:,condition_code(:,1,2))), diag(fmb_analysis.hand_scored_events(:,condition_code(:,2,2))));
    [h p] = ttest(a, b, 0.05, 'both');
    fmb_analysis.ttest.hand_scored_events = p;
    [h p] = ttest2(a, b, 0.05, 'both');
    fmb_analysis.ttest2.hand_scored_events = p;
end
if isfield(fmb_analysis, 'hand_scored_sniffs')
    a = cat(1,diag(fmb_analysis.hand_scored_sniffs(:,condition_code(:,1,1))), diag(fmb_analysis.hand_scored_sniffs(:,condition_code(:,2,1))));
    b = cat(1,diag(fmb_analysis.hand_scored_sniffs(:,condition_code(:,1,2))), diag(fmb_analysis.hand_scored_sniffs(:,condition_code(:,2,2))));
    [h p] = ttest(a, b, 0.05, 'both');
    fmb_analysis.ttest.hand_scored_sniffs = p;
    [h p] = ttest2(a, b, 0.05, 'both');
    fmb_analysis.ttest2.hand_scored_sniffs = p;
end
if isfield(fmb_analysis, 'hand_scored_false_alarms')
    a = cat(1,diag(fmb_analysis.hand_scored_false_alarms(:,condition_code(:,1,1))), diag(fmb_analysis.hand_scored_false_alarms(:,condition_code(:,2,1))));
    b = cat(1,diag(fmb_analysis.hand_scored_false_alarms(:,condition_code(:,1,2))), diag(fmb_analysis.hand_scored_false_alarms(:,condition_code(:,2,2))));
    [h p] = ttest(a, b, 0.05, 'both');
    fmb_analysis.ttest.hand_scored_false_alarms = p;
    [h p] = ttest2(a, b, 0.05, 'both');
    fmb_analysis.ttest2.hand_scored_false_alarms = p;
end
if isfield(fmb_analysis, 'hand_scored_event_failures')
    a = cat(1,diag(fmb_analysis.hand_scored_event_failures(:,condition_code(:,1,1))), diag(fmb_analysis.hand_scored_event_failures(:,condition_code(:,2,1))));
    b = cat(1,diag(fmb_analysis.hand_scored_event_failures(:,condition_code(:,1,2))), diag(fmb_analysis.hand_scored_event_failures(:,condition_code(:,2,2))));
    [h p] = ttest(a, b, 0.05, 'both');
    fmb_analysis.ttest.hand_scored_event_failures = p;
    [h p] = ttest2(a, b, 0.05, 'both');
    fmb_analysis.ttest2.hand_scored_event_failures = p;
end

fmb_analysis.version = fmba_version;

save('fmb_analysis.mat', 'fmb_analysis');
fprintf('.fmbn files in %s analyzed.\n  n files: %d.\n  analyze_fmb version %.1f\n', directory, size(filenames,2), fmba_version);
