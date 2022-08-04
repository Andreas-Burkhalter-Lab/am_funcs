%%%% removeBadFrames: eliminate stim trials during which bad 2p images were
%%%% acquired
%
% [res_out stimpar_sets_out] = removeBadFrames(res_in, stimpar_sets_in, scope_events, scopetiming, stimdur_scans, baseline_window_scans)
%
%%%% requires an excel .xlsx file with the same name as triggertiming file with '_badframes.xlsx' appended, 
%%%%        to be present in the working directory....
%%%%        ... this file should contain a first column headed by 'startframe'
%%%%        and second column headed by 'lastframe'
%%%%%%%%%%%% all frame numbers are indices within THIS PLANE,
%%%%%%%%%%%% not indices within the .sbx file or larger .tif file containing all of the
%%%%%%%%%%%% planes; the same frame numbers are deleted from each plane
% updated 2020/5/05 on thermaltake

function [res_out, stimpar_sets_out] = removeBadFrames(res_in, stimpar_sets_in, scope_events, scopetiming, stimdur_scans, baseline_window_scans)

stimpar_sets_out = stimpar_sets_in;
res_out = res_in;

bad_frames_filename = [getfname(scopetiming.triggertiming_filename), '_badframes.xlsx'];
if exist(bad_frames_filename,'file') || isfield(scopetiming,'badframes')
    ntrials = height(stimpar_sets_in);
    fprintf(['Eliminating trials containing bad frames found in ' bad_frames_filename '...\n'])
    [framelims heads] = xlsread(bad_frames_filename);
    badframes_table = table(framelims(:,1),framelims(:,2),'VariableNames',heads);
    badframes = [];
    for rr = 1:height(badframes_table)
        badframes = [badframes, badframes_table.startframe(rr):badframes_table.lastframe(rr)];
    end
    deletetrial = false(ntrials,1);
    for itrial = 1:ntrials % check if each trial (including preceding baseline period) occurred during bad frames
        timepoint_start = stimpar_sets_in.stim_onset(itrial) - baseline_window_scans; 
        timepoint_end = stimpar_sets_in.stim_onset(itrial)+stimdur_scans;
        sevents_in_range = find(scope_events.onset>timepoint_start & scope_events.onset<timepoint_end);
        sevents_in_range = [min(sevents_in_range)-1; sevents_in_range]; % include sevents that were in progress when stim started
        badframes_in_range = intersect(sevents_in_range,badframes);
        if ~isempty(badframes_in_range) % if bad frames occurred during this trial
            deletetrial(itrial) = true;
        end
        clear badframes_in_range sevents_in_range 
    end
    stimpar_sets_out(deletetrial,:) = [];
    fprintf(['Deleted ' num2str(length(find(deletetrial))) ' bad trials out of ' num2str(ntrials) ' total trials.\n'])
    res_out.removed_trials = find(deletetrial);
    res_out.badframes_table = badframes_table;
    res_out.badframes = badframes;
    res_out.bad_frames_filename = bad_frames_filename;
end