%%% get start and stop points of microscope recording within trigger timing file for imaging,   
%%%     save into file
% updated 2019/04/18 on thermaltake

function scopetiming = getScopeStartStop(triggertiming_file, savename_suffix)

plot_stimchan_events = 1; % plot stim events along with scope events, for reference; only 4th floor scope for now
fourth_floor_scope_chan_name = 'Input3'; % name of input channel for scope event timing
fourth_floor_stim_chan_name = 'Input2'; % name of input channel for stim event timing
han_lab_scope_chan_name = 'ImageSync2'; % han lab - name of input channel for scope event timing

if ~exist('triggertiming_file','var') || isempty(triggertiming_file)
    if ~isempty(find_files_of_type('abf')) % if there's abf file in this directory
        triggertiming_file = uigetfile('*.abf','Select trigger timing file.');
    else
        triggertiming_file = uigetfile('*trigdata.mat','Select trigger timing file.'); %% parsed prairie file
    end
end
if exist('savename_suffix','var')
    savename_suffix = ['_', savename_suffix];
else
    savename_suffix = '';
end

file_ext = triggertiming_file(end-3:end); 
switch file_ext
    case '.abf' % if we're trying to load an abf file
        try
            [abf_timepoints abf_sampling_interval_us abfinfo] = abfload(triggertiming_file); % abf_sampling_interval_us = sampling interval microseconds
        catch abf_mexception   %%%% if data is in .mat format instead of original .abf format
            fprintf('ABF file not readable as .abf; loading as .mat file instead.\n')
            load(triggertiming_file,'-mat')
        end
        scopetiming.scope_chan = find(strcmp(abfinfo.recChNames,han_lab_scope_chan_name));
        scopetiming.scope_chan_name = han_lab_scope_chan_name;
        timepoints_to_plot = abf_timepoints(:,scopetiming.scope_chan); % scope channel timepoints
    otherwise % 4th floor scope, prairie view
        trigdata = load_prairie_triggertiming(triggertiming_file);
        timepoints_to_plot = trigdata{:,fourth_floor_scope_chan_name}; % scope channel timepoints
        if plot_stimchan_events
            timepoints_to_plot = [timepoints_to_plot, trigdata{:,fourth_floor_stim_chan_name}]; %%% add on stim channel
        end
        scopetiming.scope_chan_name = fourth_floor_scope_chan_name;
end



plot(timepoints_to_plot)

fprintf('Zoom in to scope start then press Enter.')
k = 0;
while k ~= 1
    k = waitforbuttonpress;
end
fprintf('Click on timepoint BEFORE scope start and AFTER unwanted extraneous scope TTLs.')
startpoint = ginput(1);
scopetiming.scope_start_timepoint = startpoint(1);

fprintf('\nZoom in to scope stop then press Enter.')
k = 0;
while k ~= 1
    k = waitforbuttonpress;
end
fprintf('Click on timepoint AFTER scope end and BEFORE unwanted extraneous scope TTLs.')
stoppoint = ginput(1);
scopetiming.scope_stop_timepoint = stoppoint(1);

scopetiming.triggertiming_filename = triggertiming_file;
save([getfname(triggertiming_file) '_scopetiming', savename_suffix],'scopetiming')