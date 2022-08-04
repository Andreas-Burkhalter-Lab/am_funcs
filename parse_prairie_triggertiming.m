%%% process voltage recording from prairie view to extract stim timing, scope timing, and treadmill data
%
% parse_prairie_triggertiming(trigger_data_file_name)
%       argument can be numeric to indicate number of files to select, or can directly input a filename, 
%           or leave blank for one file and manually select
%     
% updated 2019-2-27

function parse_prairie_triggertiming(trigger_data_file_name)

treadmill_chan_1_name = 'Input0'; % first channel where treadmill voltage was collected
treadmill_chan_2_name = 'Input1'; % second channel where treadmill voltage was collected
prairie_stim_chan_name = 'Input2'; % name of input channel for visual stim timing
scope_chan_name = 'Input3'; % name of input channel for scope frame acquisition timing
filename_suffix = '_parsed';

if ~exist('trigger_data_file_name','var') || isempty(trigger_data_file_name) % select 1 file
    [trigfilename, filepath] = uigetfile('.csv','Select trigger_data_file');
    trigger_data_file_name  = fullfile(filepath, trigfilename);
    nfiles = 1;
    savename = [getfname(trigger_data_file_name), filename_suffix];
    trigger_data_file_name = {trigger_data_file_name}; 
elseif isnumeric(trigger_data_file_name) % multiple files
    nfiles = trigger_data_file_name;
    clear trigger_data_file_name
    for ifile = 1:nfiles
         [trigfilename, filepath] = uigetfile('.csv',['Select trigger_data_file ', num2str(ifile)]);
         trigger_data_file_name{ifile,1}  = fullfile(filepath, trigfilename); 
    end
    pathparts = strsplit(pwd,'\');
    savename = [pathparts{end-1} '_'  pathparts{end} '_trigdata'];
else % filename specified in argument
    nfiles = 1;
    savename = [getfname(trigger_data_file_name), filename_suffix];
    trigger_data_file_name = {trigger_data_file_name}; 
end

trigdata = table;
for ifile = 1:nfiles
    this_filename = trigger_data_file_name{ifile};
    trigdata_thisfile = load_prairie_triggertiming(this_filename);
    % get treadmill data if available
    if all(ismember({treadmill_chan_1_name, treadmill_chan_2_name}, trigdata_thisfile.Properties.VariableNames)) %% if treadmill channels were recorded
        speedtable = parse_treadmill_velocity([trigdata_thisfile{:,treadmill_chan_1_name}, trigdata_thisfile{:,treadmill_chan_2_name}]);
        trigdata_thisfile.locm_forw_mps = speedtable(:,5); % only save the smoothed velocity in meters per second
        trigdata_thisfile(:,treadmill_chan_1_name) = [];
        trigdata_thisfile(:,treadmill_chan_2_name) = [];
    end
    trigdata_thisfile.Time_ms = [];
    trigdata_thisfile{:,prairie_stim_chan_name} = uint8(trigdata_thisfile{:,prairie_stim_chan_name}); % decrease file size
    trigdata_thisfile{:,scope_chan_name} = uint8(trigdata_thisfile{:,scope_chan_name}); % decrease file size
    
    trigdata = [trigdata; trigdata_thisfile]; % concatenate data from this file to the full table
end

save(savename,'trigdata','trigger_data_file_name')
fprintf(['Trigger timing data parsed and saved into ' savename '.mat \n'])