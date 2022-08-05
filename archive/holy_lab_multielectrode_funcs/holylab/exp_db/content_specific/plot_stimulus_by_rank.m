function out_handle = plot_stimulus_by_rank(analysis_file_list, stimulus, varargin)
% function plot_stimulus_by_rank plots a grid-style summary of responses to
% a set of stimuli ranked based on their responses to a particular stimulus
% Required input:
%     analysis_file_list: A full-path cell array of strings of the
%                         filenames of 'analysis.mat' to be analyzed
%     stimulus:  can be:  (a) a verbatim character matrix containing the name
%                             of the stimulus for rank
%                         (b) an integer in the range of used valves which
%                             is the valvenumber to be analyzed
% Optional input: (in progress)
%   option structure containing one or more of these fields:
%     'axis' = axis handle for plotting, otherwise data will be plotted in
%     an axis in a new figure
% Output:
%     out_handle: handle of the axis to which data were plotted
%

% Copyright 2008 Julian P. Meeks (Timothy Holy laboratory)
% Revision history:
%       2008_02_08: Wrote it (JPM)
%

%% Section 1: Input checking
% first variable (analysis_file_list)-------------------------------------
if ~iscell(analysis_file_list)
    if ~ischar(analysis_file_list)
        errormsg = 'analysis_file_list must be a string or cell array of strings';
        error(errormsg);
    else
        analysis_file_list = {analysis_file_list}; % if single string input, cell-it
    end
else
    for file_idx = 1:size(analysis_file_list,2)
        if ~ischar(analysis_file_list{file_idx})   % fail if any of cells don't have a string
            errormsg = 'analysis_file_list must be a cell array of strings and strings only';
            error(errormsg);
        end
        if ~exist(analysis_file_list{file_idx}, 'file')  % fail if any files don't exist
            errormsg = [analysis_file_list{file_idx} ' is not a valid filename'];
            error(errormsg);
        end
    end
end
%  ----- end checking for 'analysis_file_list' variable ------------------

% second variable (stimulus)----------------------------------------------
if iscell(stimulus)
    if size(stimulus,2)~=1         % only 1 stimulus okay right now
        errormsg = 'only one stimulus may be plotted at a time';
        error(errormsg);
    else
        stimulus = stimulus{1};     % un-cell any cell-ed input
    end
else
    if ~isa(stimulus, 'numeric') && ~ischar(stimulus) % check for proper format
        errormsg = 'stimulus must be specified by unique string tag or valve number';
        error(errormsg);
    end
end
% ----- end checking for 'stimulus' variable -----------------------------

% varargin checking ------------------------------------------------------
 % step 1 (currently the only step) = Check for input axis
if nargin > 2
    options = varargin{1};
    if ~isstruct(options)        % check to make sure is in structure format
        errormsg = 'varargin must be an options structure, consult help';
        error(errormsg);
    end
    if isfield(options, 'axis')
        plot_axis = options.axis;
    else
        plot_axis = 0;
    end
end

%% Establish plot axis
if ~exist('plot_axis', 'var');
    plot_axis = 0;
end
if plot_axis == 0;  % if input axis not present, make a new figure and plot axis
    mainfig = figure;
    plot_axis = axes('parent', mainfig, 'position', [0.1 0.1 0.9 0.9], 'units', 'normalized');
end
%% Load and check 'analysis.mat' files
% load data files
data = gather_analysis(analysis_file_list);   % returns structure array with primary contents
                                              %  in structure called 'analysis'
                                              %  (i.e. data{1}.analysis contains the values for
                                              %  the first data point)
                                              
% convert stimulus to an array of the proper valves to check
if ~isa(stimulus,'numeric')        % if a valve label is given, find the appropriate valves
    for data_idx = 1:size(data,2)      % NOTE: size of data array may be larger than size of file list
        matches = strmatch(stimulus, data{data_idx}.analysis.valvelabels);
        if ~isempty(matches)  % if the string exists in the valvelabels list
            valve_indices(data_idx) = matches;
        else
            valve_indices(data_idx) = NaN;
        end
    end
else
    for data_idx = 1:size(data,2)  % check to make sure there are assigned valves of a
        active_valves = 0;         % valve input number
        for valve_idx = 1:size(data{data_idx}.analysis.valvelabels,2)
            active_valves(valve_idx) = ischar(data{data_idx}.analysis.valvelabels{valve_idx});
        end
        n_valves(data_idx) = sum(active_valves);
    end
    if sum(n_valves) == 0;
        errormsg = 'none of the experiments used the valve number supplied';
        error(errormsg);
    end
    valve_indices = zeros(size(data));  % if a valve number was given, each file will have
    valve_indices(:) = stimulus;           % that valve checked (regardless of its label)
end

if sum(valve_indices) == 0 || isnan(sum(valve_indices))   % throw an error if this didn't produce any valves
    errormsg = 'no valves with that label were found in the data list, please check spelling';
    error(errormsg);
end

%% Find cells which were responsive (statistical p value < 0.06) to the stimulus
for data_idx = 1:size(data,2)
    if data{data_idx}.analysis.t_test_pvalue(valve_indices(data_idx)) < 0.06
        if ~exist('responsive_data', 'var')
            responsive_data = data{data_idx};
        else
            responsive_data(end+1) = data{data_idx};
        end
    end
end

%% Rank these cells based on their firing rate to the key stimulus
ranked_data = cell(size(responsive_data));
delta_r_orig = zeros(1, size(responsive_data,2));
for resp_idx = 1:size(responsive_data,2)  % make a simple array of the mean delta_r vals at play
    delta_r_orig(resp_idx) = responsive_data(resp_idx).analysis.delta_r_mean(stimulus);
end
[sorted_delta_r, sorted_indices] = sort(delta_r_orig,2,'descend');  % sort them
for resp_idx = 1:size(responsive_data,2)              % use sorted_indicies to re-rank data
    ranked_data(resp_idx) = {responsive_data(sorted_indices(resp_idx))};
end

%% Plot the data
% set up the sub-axes (left column = reference image, right block = normalized data grid)
sub_axes = SplitHoriz([1/max(n_valves) 2/max(n_valves)],[1 0 1], plot_axis);
reference_plot = sub_axes(1);
normalized_plot = sub_axes(2);

% Plot reference_plot
reference_img(:,1) = sorted_delta_r(:);
imagesc(reference_img,'parent', reference_plot);
set(reference_plot, 'clim', [-mean(reference_img(:,1))-std(reference_img(:,1)),...
                              mean(reference_img(:,1))+std(reference_img(:,1))]);

% Set up normalized image
for ranked_idx = 1:size(ranked_data,2)
    if mean(ranked_data{ranked_idx}.analysis.delta_r_mean) > 0
        norm_img(ranked_idx,:) = ranked_data{ranked_idx}.analysis.delta_r_mean./...
                                 mean(ranked_data{ranked_idx}.analysis.delta_r_mean);
    elseif mean(ranked_data{ranked_idx}.analysis.delta_r_mean) == 0
        norm_img(ranked_idx,:) = ranked_data{ranked_idx}.analysis.delta_r_mean./...
                                 std(ranked_data{ranked_idx}.analysis.delta_r_mean);
    elseif mean(ranked_data{ranked_idx}.analysis.delta_r_mean) < 0
        norm_img(ranked_idx,:) = ranked_data{ranked_idx}.analysis.delta_r_mean./...
                                 (-1*mean(ranked_data{ranked_idx}.analysis.delta_r_mean));
    end
end
imagesc(norm_img, 'parent', normalized_plot);
set(normalized_plot, 'clim', [-10 10]);

%% Assign output
out_handle = plot_axis;

end