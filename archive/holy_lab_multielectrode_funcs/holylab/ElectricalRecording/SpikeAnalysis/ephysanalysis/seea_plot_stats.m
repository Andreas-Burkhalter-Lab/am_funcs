function seea_plot_stats(analysis_in, test, category)
% SEEA_PLOT_STATS creates a 3-d plot of p values in an analysis struct
% Syntax: seea_plot_stats(analysis_in, test, category)
%         each field is required
%  - 'analysis_in' can be a structure variable or directory containing a
%       *analysis.mat file
%  - 'test' must be a character array containing the name of the desired
%       statistical test used to create the p values
%  - 'category' must be a character array including the name of the field
%       you wish to view ('delta_r', 'fr', 'base_r', etc.)
%
% Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)

% Revision History:
% 2008_05_09 Wrote it (JPM)

%% Step 1: Check inputs
%  1a: if directory supplied, check to make sure we have data to plot, if
%      so assign data to analysis_in variable.  If not directory, check to
%      make sure that analysis_in is a structure variable;
if ischar(analysis_in)
    if ~exist(analysis_in, 'dir')
        errormsg = '''analysis_in'' must be a directory name containing a *analysis.mat file.';
        error(errormsg);
    elseif size(dirbyname([analysis_in filesep '*analysis.mat']),2)>1
        errormsg = [analysis_in 'directory contains more than one *analysis.mat file.'];
        error(errormsg);
    elseif size(dirbyname([analysis_in filesep '*analysis.mat']),2)<1
        errormsg = [analysis_in 'directory does not contain an *analysis.mat file.'];
        error(errormsg);
    else
        filename = dirbyname([analysis_in filesep '*analysis.mat']);
        load([analysis_in filesep filename{1}]);
        if ~exist('analysis', 'var')
            errormsg = 'The input structure in *analysis.mat file should be named ''analysis''.';
            error(errormsg);
        else
            analysis_in = analysis;
        end
    end
end
if ~isstruct(analysis_in)
    errormsg = '''analysis_in'' must be a string containing the name of a directory or an ''analysis'' structure variable';
    error(errormsg);
end

% % 1b: Check to make sure the desired 'category' is present in analysis struct
% if ~isfield(analysis_in, category)
%     errormsg = sprintf('The %s category does not exist in the analysis_in structure.',category);
%     error(errormsg);
% end

% 1b: Check to make sure the desired 'test' is preset in the analysis struct
if ~isfield(analysis_in, [category '_' test])
    errormsg = ['The ' test ' test does not exist for ' category ' in the analysis_in structure.'];
    error(errormsg);

    % 1c: If certain categories are chosen, we have to take into consideration
    %     that the data represent comparisons between certain valves (saved
    %     into analysis substructure; for example:
    %     fem_selectivity_f_valves/m_valves
else
    switch [category '_' test]
        case 'fem_selectivity_ttest'
            valvenames = analysis_in.valvelabels(analysis_in.fem_selectivity_f_valves); % sloppy for now
            yticks = 1:size(valvenames,2);
    end
    if ~exist('valvenames', 'var')
        valvenames = analysis_in.valvelabels;
        yticks = 1:size(valvenames,2);
    end
end


%% Step 2 Load in the desired data set
fieldtoplot = [category '_' test];
field_names = fieldnames(analysis_in.(fieldtoplot));
fieldname_size = size(field_names,1);
valves = analysis_in.valvelabels;
n_valves_toplot= size(analysis_in.(fieldtoplot).(field_names{1}),1);
n_cells = size(analysis_in.(fieldtoplot).(field_names{1}),2);

% Initialize variable space [valves x field_names x cells]
toplot = ones(n_valves_toplot, fieldname_size, n_cells);

for c_idx = 1: n_cells
    for f_idx = 1:fieldname_size
        toplot(:,f_idx,c_idx) = 1./(analysis_in.(fieldtoplot).(field_names{f_idx})(:,c_idx));
    end
    hfig(c_idx) = figure('units','pixels','position',[100 100 900 700], 'name', 'seea_plot_stats');
    hax(c_idx) = imagesc(log10(toplot(:,:,c_idx)));  % plotting the inverse log10 of the p value
    if isfield(analysis_in, 'defined_subintervals')
        if iscell(analysis_in.defined_subintervals)
            xticks = 1:size(analysis_in.defined_subintervals,2);
            xticklabels = analysis_in.defined_subintervals;
        elseif analysis.defined_subintervals == 'all'
            xticks = 1:10:size(field_names,2);
            xticklabels = field_names(1:10:size(field_names,2));
        else
            xticks = 1;
            xticklabels = {analysis_in.defined_subintervals};
        end
    end
    set(get(hax(c_idx),'parent'),'title', text('string',analysis_in.id),...
        'ytick', yticks,...
        'yticklabel', valvenames,...
        'xtick', xticks,...
        'xticklabel', xticklabels);
    
    chax(c_idx) = colorbar;
end



end
