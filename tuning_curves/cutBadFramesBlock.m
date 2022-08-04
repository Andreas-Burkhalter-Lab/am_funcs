%%%% cutBadFramesBlock: cut out from scope_events all frames that occurred before a bad frame
%%%% listed in the badframes file
%
% [res_out, scope_events_out] = cutBadFramesBlock(res_in, scope_events_in, scopetiming)
%
%%%%%% discard whichever section is smaller: the block before the last bad frame 
%%%%%%      or the block after the first bad frame
%%%% requires an excel .xlsx file with the same name as scopeteiming_file, but with 'scopetiming'
%%%% replaced with 'badframes', will need to be present in the working
%%%% dir.... this file should contain a first column headed by 'startframe'
%%%% and second columnd headed by 'endframe'
%%%%%%%%%%%% all frame numbers are indices within THIS PLANE,
%%%%%%%%%%%% not indices within the .sbx file containing all of the
%%%%%%%%%%%% planes; the same frame numbers are deleted from each plane
% updated 5/17/18 by Andrew Meier on thermaltake

function [res_out, scope_events_out] = cutBadFramesBlock(res_in, scope_events_in, scopetiming)

scope_events_out = scope_events_in;
res_out = res_in;

bad_frames_filename = [fileparts(scopetiming.abf_filename), filesep, getfname(scopetiming.abf_filename), '_badframes.xlsx'];
if exist(bad_frames_filename,'file') || isfield(scopetiming,'badframes')
    fprintf(['Eliminating all scope events for this plane occurring before any bad frame contained in ' bad_frames_filename '...\n'])
    [framelims heads] = xlsread(bad_frames_filename);
    badframes_table = table(framelims(:,1),framelims(:,2),'VariableNames',heads);
    badframes = [];
    for rr = 1:height(badframes_table)
        badframes = [badframes, badframes_table.startframe(rr):badframes_table.lastframe(rr)];
    end

    %%% determine whether frames before or after the bad frames should be kept
    nframesraw = height(scope_events_in);
    firstbadframe = min(badframes);
    lastbadframe = max(badframes);
    if firstbadframe > nframesraw - lastbadframe % if the period before bad frames is longer than the period after bad frames
        res_out.frame_range_removed = [firstbadframe nframesraw];
        res_out.frame_range_kept = [1 firstbadframe-1];
        scope_events_out = scope_events_out(1:firstbadframe-1,:); % delete frames occuring during or after the first bad frame
    else % if the period after bad frames is longer than the period before bad frames
        res_out.frame_range_removed = [1 lastbadframe];
        res_out.frame_range_kept = [lastbadframe+1 nframesraw];
        scope_events_out = scope_events_out(lastbadframe+1:end,:); % delete frames occuring before or during a bad frame
    end
    res_out.n_frames_removed = res_out.frame_range_removed(2) - res_out.frame_range_removed(1) + 1;

    fprintf(['Deleted ' num2str(res_out.n_frames_removed) ' frames ([' num2str(res_out.frame_range_removed) ']) out of ' num2str(nframesraw) ' total frames.\n'])
    
    res_out.nframesraw = nframesraw;
    res_out.nframes = nframesraw - res_out.n_frames_removed;
    res_out.badframes_table = badframes_table;
    res_out.badframes = badframes;
    res_out.bad_frames_filename = bad_frames_filename;
end