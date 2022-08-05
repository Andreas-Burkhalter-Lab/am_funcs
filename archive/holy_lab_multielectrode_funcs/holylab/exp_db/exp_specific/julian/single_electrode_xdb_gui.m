function entry = single_electrode_xdb_gui(xdb_id) 
%SINGLE_ELECTRODE_XDB_GUI is a gui-based function which will allow a user
%to view and edit information from single electrode .merec recordings and
%subsequent snippets and spike sorting data.
%
% single_electrode_xdb_gui takes a single argument xdb_id which is a string
%       containing the xdb_id, without full path, and without the .xdb extension.
% The user may supply xdb_id 'new' in order to create a template for a new
%       .xdb entry, which he/she will have the option to save after
%       supplying a minimum (required) amount of information
%
% Copyright 2007 Julian P. Meeks (Timothy Holy Laboratory)
%
%   
%
% See also SNIPOPTIONS, SAVE_XDB_ENTRY, SAVE_XDB_COOKIE

%% Load in previously-saved .xdb data
entry = struct;
if strmatch(xdb_id, 'new', 'exact') % if the user wants a new, blank .xdb entry
    xdb_id = save_xdb_entry(entry); % create new xdb_id with empty entry
end

username = xdb_id(1:strfind(xdb_id, '_')-1); % extract username from xdb_id
entry = load(strcat(['/usr/lab/exp_db/' username filesep xdb_id '.xdb']), '-mat'); % load entry

%% Load in previously-saved .analysis data

% yet to be implemented.  A .analysis file will presumably be the "data"
% portion of Rebecca's dcell structure (i.e. spike times, sorting links,
% other analysis linked to intervals, etc).  This should be the portion of
% the dcell structure that has the potential to balloon into larger sizes.


%% Initialize GUI panels and fields

% Load a figure 'mainfig'
mainfig = figure(1);
set(mainfig, 'Units', 'normalized', 'Name', 'Single Electrode .xdb GUI',...
        'Position', [0 0 1 1]);

% Initialize version_panel (simple: contains only the version number, which
% should be hard-coded at initiation and updated when .xdb receives updates
version_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.01, 0.92, 0.18,0.07],...
                        'BackgroundColor', [.8 .34 .45]);
version_num_text = uicontrol('parent', version_panel, 'style', 'text', 'string', '.XDB Version Number:',...
                        'Units', 'normalized', 'position', [0.01, 0.21, 0.78, 0.68], 'FontSize', 14,...
                        'BackgroundColor', [.8 .34, .45],'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'left');
set(version_num_text, 'FontUnits', 'normalized');
if ~isfield(entry, 'version'); entry.version = 1.0; end     % set manually if not already present
version_num = uicontrol('parent', version_panel, 'Units', 'normalized', 'Position', [0.81, 0.21, 0.18, 0.68],...
                        'FontSize', 14, 'BackgroundColor', [.8 .34 .45], 'FontWeight', 'bold',...
                        'String', num2str(entry.version), 'style', 'text');
% ----------/version_panel------------------------------------------------

% Initialize investigator_panel ------------------------------------------
    % ASSERT: username must be one of the investigators
    if ~isfield(entry, 'investigator')
        entry.investigator = {username};
    end 
    %--------------------------------------------------
investigator_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.01, 0.82, 0.18,0.07],...
                        'BackgroundColor', [.23 .74, .45]);
investigator_text = uicontrol('parent', investigator_panel, 'style', 'text', 'string', 'Investigator(s):',...
                        'Units', 'normalized', 'position', [0.01, 0.21, 0.78, 0.68], 'FontSize', 14,...
                        'BackgroundColor', [.23 .74, .45],'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'left');
set(investigator_text, 'FontUnits', 'normalized');
investigator_names = uicontrol('parent', investigator_panel, 'Units', 'normalized', 'Position', [0.51, 0.31, 0.48, 0.68],...
                        'FontSize', 14, 'BackgroundColor', [.23 .74, .45], 'FontWeight', 'bold',...
                        'String', entry.investigator, 'style', 'popupmenu');
investigator_edit_btn = uicontrol('parent', investigator_panel, 'Units', 'normalized', 'Position', [0.01, 0.01, 0.78, 0.38],...
                        'FontSize', 14,...
                        'String', 'Add/Edit investigator(s)', 'style', 'pushbutton', 'callback', @edit_investigator);
% ----------/investigator_panel--------------------------------------------

% Initialize start_date_panel ---------------------------------------------
%   Purpose: show/edit the start date of the experiment
start_date_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.01, 0.72, 0.18,0.07],...
                        'BackgroundColor', [.83 .44, .15]);
start_date_text1 = uicontrol('parent', start_date_panel, 'style', 'text', 'string', 'Start date:',...
                        'Units', 'normalized', 'position', [0.01, 0.21, 0.28, 0.68], 'FontSize', 14,...
                        'BackgroundColor', [.83, .44, .15], 'fontweight', 'bold',...
                        'HorizontalAlignment', 'left');
start_date_text2 = uicontrol('parent', start_date_panel, 'style', 'text', 'string', 'YEAR/MONTH/DAY',...
                        'Units', 'normalized', 'position', [0.31, 0.21, 0.68, 0.68], 'FontSize', 14,...
                        'BackgroundColor', [.83, .44, .15], 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'left');                  
% Initialize the start date field ----------------------------------------
if ~isfield(entry, 'start_date')
    entry.start_date = [str2num(datestr(now, 'yyyy'))...
        str2num(datestr(now, 'mm'))...
        str2num(datestr(now, 'dd'))];
end
 setappdata(start_date_panel, 'start_date_month', entry.start_date(2));
 setappdata(start_date_panel, 'start_date_day', entry.start_date(3));
 setappdata(start_date_panel, 'start_date_year', entry.start_date(1));
% --start date fields are set --------------------------------------------
start_date_year = uicontrol('parent', start_date_panel, 'style', 'edit',...
                              'string', getappdata(start_date_panel, 'start_date_year'),...
                              'units', 'normalized', 'fontsize', 12, 'position', [.27 .01 .25 .48],...
                              'callback', @update_start_year);
start_date_month = uicontrol('parent', start_date_panel, 'style', 'edit',...
                              'string', getappdata(start_date_panel, 'start_date_month'),...
                              'units', 'normalized', 'fontsize', 12, 'position', [.57 .01 .15 .48],...
                              'callback', @update_start_month);
start_date_day = uicontrol('parent', start_date_panel, 'style', 'edit',...
                              'string', getappdata(start_date_panel, 'start_date_day'),...
                              'units', 'normalized', 'fontsize', 12, 'position', [.82 .01 .15 .48],...
                              'callback', @update_start_day);
%----------------------/start_date_panel----------------------------------

% Initialize data_locations_panel
%   This panel will allow the user to browse for specific files to add to
%   the data_locations field
data_locations_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.01, 0.60, 0.18,0.1],...
                        'BackgroundColor', [.83 .44, .95]);
data_locations_text1 = uicontrol('parent', data_locations_panel, 'style', 'text', 'string', 'Data File(s):',...
                        'Units', 'normalized', 'position', [0.01, 0.41, 0.33, 0.58], 'FontSize', 12,...
                        'BackgroundColor', [.83 .44, .95],'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'left');
if isfield(entry, 'data_locations')
    if ~isempty(entry.data_locations)
        for idx = 1:size(entry.data_locations,2)
            data_locations(idx) = {strip_path(entry.data_locations{idx})};
        end
        data_locations_popup = uicontrol('parent', data_locations_panel, 'Units', 'normalized', 'Position', [0.01, 0.01, 0.98, 0.48],...
                        'FontSize', 10, 'BackgroundColor', [.83 .44, .95], 'FontWeight', 'normal',...
                        'String', data_locations,...
                        'style', 'popupmenu');
    end
else
    data_locations_popup = uicontrol('parent', data_locations_panel, 'Units', 'normalized', 'Position', [0.01, 0.01, 0.98, 0.48],...
                        'FontSize', 10, 'BackgroundColor', [.83 .44, .95], 'FontWeight', 'normal',...
                        'String', ' ', 'style', 'popupmenu');    
end
data_locations_addbtn = uicontrol('parent', data_locations_panel, 'Units', 'normalized', 'Position', [0.36, 0.51, 0.31, 0.48],...
                        'FontSize', 12,'String', 'Add files(s)', 'style', 'pushbutton',...
                        'callback', @add_data_locations);
data_locations_removebtn = uicontrol('parent', data_locations_panel, 'Units', 'normalized', 'Position', [0.68, 0.51, 0.31, 0.48],...
                        'FontSize', 12,'String', 'Remove file', 'style', 'pushbutton',...
                        'callback', @remove_data_locations);
% ---------------------/data_location_panel-------------------------------                    

% initialize tags_panel to allow user to select the specific tags associated
%   with this particular experiment
% tags.mat is required (@ /matlabfunc/exp_db/exp_specific/julian/tags.mat)
%   --> this file contains the allowed tags to be put into the
%       tags_base_list

% load base tags
    tags = cell(1);
    load /mnt/julian_004/julian/matlabfunc/exp_db/exp_specific/julian/tags.mat
    if isfield(entry, 'tags')
        if ~isempty(entry.tags)
            for idx_tags = 1:size(entry.tags,2)
                matched = strmatch(entry.tags{idx_tags}, tags);
                if ~isempty(matched)
                    tags(matched) = [];
                end
            end
        end
    end
% ASSERT: <tags> contains the list of available tagnames in a cell array of strings

tags_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.01, 0.11, 0.18,0.48],...
                        'BackgroundColor', [.83 .94, .05]);
tags_text1 = uicontrol('parent', tags_panel, 'style', 'text', 'string', 'Tags:',...
                        'Units', 'normalized', 'position', [0.01, 0.93, 0.33, 0.06], 'FontSize', 12,...
                        'BackgroundColor', [.83 .94, .05],'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'left');
tags_text2 = uicontrol('parent', tags_panel, 'style', 'text', 'string', 'Available:',...
                        'Units', 'normalized', 'position', [0.11, 0.88, 0.33, 0.06], 'FontSize', 12,...
                        'BackgroundColor', [.83 .94, .05],'FontWeight', 'normal', ...
                        'HorizontalAlignment', 'left');
tags_text3 = uicontrol('parent', tags_panel, 'style', 'text', 'string', 'Selected:',...
                        'Units', 'normalized', 'position', [0.61, 0.88, 0.33, 0.06], 'FontSize', 12,...
                        'BackgroundColor', [.83 .94, .05],'FontWeight', 'normal', ...
                        'HorizontalAlignment', 'left');
tags_base_list = uicontrol('parent', tags_panel, 'style', 'listbox', 'string', tags,...
                        'Units', 'normalized', 'position', [0.01, 0.01, 0.43, 0.89], 'FontSize', 10,...
                        'BackgroundColor', [1 1, 1],'FontWeight', 'normal', 'max', size(tags,2));
if ~isfield(entry, 'tags')
    entry.tags = [];
end
if isempty(entry.tags)
    tags_selected_list =  uicontrol('parent', tags_panel, 'style', 'listbox', 'string', 'empty',...
                        'Units', 'normalized', 'position', [0.56, 0.01, 0.43, 0.89], 'FontSize', 10,...
                        'BackgroundColor', [1 1, 1],'FontWeight', 'normal', 'max', size(tags,2));
else
    tags_selected_list =  uicontrol('parent', tags_panel, 'style', 'listbox', 'string', entry.tags,...
                        'Units', 'normalized', 'position', [0.56, 0.01, 0.43, 0.89], 'FontSize', 10,...
                        'BackgroundColor', [1 1, 1],'FontWeight', 'normal', 'max', size(tags,2));
end                 
add_tag_btn = uicontrol('parent', tags_panel, 'Units', 'normalized', 'Position', [0.45, 0.61, 0.1, 0.18],...
                        'FontSize', 10,'String', '->', 'style', 'pushbutton',...
                        'callback', {@add_to_selected,tags_base_list, tags_selected_list, 'tags'});
remove_tag_btn = uicontrol('parent', tags_panel, 'Units', 'normalized', 'Position', [0.45, 0.31, 0.1, 0.18],...
                        'FontSize', 10,'String', '<-', 'style', 'pushbutton',...
                        'callback', {@remove_from_selected,tags_base_list, tags_selected_list, 'tags'});
% -----------------/tags_panel--------------------------------------------

% Initialize save_panel to save the current data to the .xdb file
save_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.01, 0.01, 0.18,0.09],...
                        'BackgroundColor', [.23 .94, .95]);
save_btn = uicontrol('parent', save_panel, 'Units', 'normalized', 'Position', [0.01, 0.01, 0.98, 0.98],...
                        'FontSize', 14,'String', 'SAVE TO .XDB', 'style', 'pushbutton',...
                        'callback', @save_entry);
% -----------------/save_panel--------------------------------------------

% Initialize tag_choice_panel: small panel to contain the popup menu for
% the middle section of the gui, which will add a specific entry form
% depending on the type of field you wish to modify.  As a default, the
% 'stimulus', 'experimental_animal', and 'analysis' options will be
% available.
tag_choice_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.21, 0.93, 0.18,0.06],...
                        'BackgroundColor', [.73 .74, .35], 'title', 'View/Edit Field Options',...
                        'fontsize', 12);
tag_choice_popup = uicontrol('parent', tag_choice_panel, 'Units', 'normalized', 'Position', [0.11, 0.11, 0.78, 0.78],...
                        'FontSize', 12, 'BackgroundColor', [.73 .74, .35], 'FontWeight', 'normal',...
                        'String', [{'stimulus'}, {'experimental_animal'}, {'analysis'} entry.tags],...
                        'style', 'popupmenu', 'callback', @display_field_options_panel);
% --------------------/tag_choice_panel
                    
%% Define callback functions

% edit_investigator BUTTON to allow user to input investigator names (4 max
% at this time)
    function result = edit_investigator(hObject, eventdata)
        if isfield(entry, 'investigator')
            names = entry.investigator;
            emptychar = cell([4 2]-size(names));
            emptychar(1:size(emptychar,1)) = {''};
            names = cat(1, names, emptychar);
            names = inputdlg({'Enter investigator name(s):', ' ', ' ', ' '}, 'Add/Edit Investigator(s)', 1, names);
        else
            names = cell([4 1]);
            names(1) = {username};
            names = inputdlg({'Enter investigator name(s):', ' ', ' ', ' '}, 'Add/Edit Investigator(s)', 1, names);
        end
        while isempty(names)
            names = inputdlg({'Enter investigator name(s):', ' ', ' ', ' '}, 'Add/Edit Investigator(s)', 1);
        end
        for idx = 1:size(names,1)
            if isempty(names{idx})
                names(idx:size(names,1)) = [];
               break
            end
        end
        entry.investigator = names;
        set(investigator_names, 'string', entry.investigator, 'horizontalalignment', 'right');
    end
%----------------/edit_investigator----------------------------------------

% update_start_year: to update just the start year using the
%   start_date_year edit text field of the gui
    function update_start_year(hObject, eventdata)
        new_str = get(hObject, 'string');
        if isnan(str2double(new_str))
            errordlg('Input must be numeric (4-digit year)', 'Error with year input');
            set(hObject, 'string', getappdata(start_date_panel, 'start_date_year'));
            return;
        elseif length(new_str)~= 4
            errordlg('Input must be a 4-digit year (i.e. 2007 not 07', 'Error with year input');
            set(hObject, 'string', getappdata(start_date_panel, 'start_date_year'));
            return;
        end

        setappdata(start_date_panel, 'start_date_year', new_str);
        entry.start_date(1) = str2num(new_str);
    end
% ------------------------/update_start_year------------------------------            

% update_start_month: to update just the start month using the
%   start_date_month edit text field of the gui
    function update_start_month(hObject, eventdata)
        new_str = get(hObject, 'string');
        if isnan(str2double(new_str))
            errordlg('Input must be numeric (2-digit month)', 'Error with month input');
            set(hObject, 'string', getappdata(start_date_panel, 'start_date_month'));
            return;
        elseif length(new_str)~= 2
            errordlg('Input must be a 2-digit month (i.e. 03 not March', 'Error with day input');
            set(hObject, 'string', getappdata(start_date_panel, 'start_date_month'));
            return;
        end

        setappdata(start_date_panel, 'start_date_month', new_str);
        entry.start_date(2) = str2num(new_str);
    end
% ------------------------/update_start_month------------------------------              

% update_start_day: to update just the start day using the
%   start_date_day edit text field of the gui
    function update_start_day(hObject, eventdata)
        new_str = get(hObject, 'string');
        if isnan(str2double(new_str))
            errordlg('Input must be numeric (2-digit day)', 'Error with day input');
            set(hObject, 'string', getappdata(start_date_panel, 'start_date_day'));
            return;
        elseif length(new_str)~= 2
            errordlg('Input must be a 2-digit day (i.e. 03 not 3', 'Error with day input');
            set(hObject, 'string', getappdata(start_date_panel, 'start_date_day'));
            return;
        end

        setappdata(start_date_panel, 'start_date_day', new_str);
        entry.start_date(3) = str2num(new_str);
    end
% ------------------------/update_start_day------------------------------

% add_data_locations: pushbutton callback to call UIGetFiles to allow user
% to select files using a GUI interface.  
% TO DO: ADD EDIT/DELETE FILES option
    function add_data_locations(hObject, eventdata)
        new_files = UIGetFiles('*.*', 'Please select the files you want to add:', pwd);
        for idx_new_files = 1:size(new_files, 2)
            if isfield(entry, 'data_locations')
                if isempty(strmatch(new_files{idx_new_files}, entry.data_locations))
                    entry.data_locations(end+1) = new_files(idx_new_files);
                end
            else
                entry.data_locations(1) = new_files(idx_new_files);
            end
        end
        for idx_datalocs = 1:size(entry.data_locations,2)
            location_nopath(idx_datalocs) = {strip_path(entry.data_locations{idx_datalocs})};
        end
            set(data_locations_popup, 'String', location_nopath, 'fontsize', 10, 'fontweight', 'normal');
    end
% ----------/add_data_locations-------------------------------------------

% remove_data_locations takes the currently selected file from the data_locations_popup
% menu and removes it from entry.data_locations (and the popup list)
    function remove_data_locations(hObject, eventdata)
        value = get(data_locations_popup, 'value');
        entry.data_locations(value) = [];
        if size(entry.data_locations,2)~=0
            for idx = 1:size(entry.data_locations,2)
                strings(idx) = {strip_path(entry.data_locations{idx})};
            end
        else
            strings = ' ';
        end
        set(data_locations_popup, 'string', strings, 'value', 1);
    end
% ------------/remove_data_locations--------------------------------------


% save_entry will take the current entry, validate and save it
    function save_entry(hObject, eventdata)
        [passed errors] = validate_xdb_entry_var(entry);
        if passed
            save_xdb_entry(entry, struct('id', entry.id));
            if isappdata(tag_choice_panel,'cur_panel_saved')
                setappdata(tag_choice_panel,'cur_panel_save',1);
            end
        else
            errordlg(errors);
        end
    end
% ----------------------/save_entry---------------------------------------

% display_field_options_panel will call up a specific panel in the center
% region of the gui to allow the user to use the gui to input the
% parameters for each experiment.
% TO DO: in cases where an xdb entry contains physiology-associated
% information which can be obtained from the file header, allow the user to
% view the info (but not necessarily edit it) using this same drop-down
% menu
    function display_field_options_panel(hObject, eventdata)
        string = get(hObject,'string');
        value = get(hObject, 'value');
        choice = string{value};
        if isappdata(tag_choice_panel, 'current_panel')
            if ~isempty(strmatch(choice,getappdata(tag_choice_panel, 'current_panel')))
                return;
            end
        end
        switch choice
            case 'stimulus'
                display_stimulus_fieldoptions(hObject, eventdata);
            case 'experimental_animal'
                display_experimental_animal_fieldoptions(hObject, eventdata);
        end
        setappdata(tag_choice_panel, 'current_panel', choice);
    end
% ----------------/display_field_options_panel----------------------------

%% DEFINE SPECIFIC GUI PANELS for XDB SUBFIELDS

% display_stimulus_fieldoptions simply opens a panel with the stimulus field options
    function display_stimulus_fieldoptions(hObject, eventdata)
        % first ensure you don't delete unsaved data currently onscreen
        if ~screen_fieldoption_panel(hObject, eventdata)==1
            return;
        end

        % initialize the stimulus structure and tag_choice_panel 'current_stimulus',etc
        if ~isfield(entry, 'stimulus')
            entry.stimulus = {struct('category', {''}, 'identity', {''}, 'concentration', {[]}, ...
                'conc_unit', {''}, 'duration', {[]}, 'user_tag', {''}, 'strain', {''}, 'sex', {''})};
            n_stimuli = 1;
            setappdata(tag_choice_panel,'n_stimuli', n_stimuli);          % store the number of stimuli in this experiment
            current_stimulus = 1;
            setappdata(tag_choice_panel,'current_stimulus', current_stimulus);   % store the current stimulus being displayed/edited
            n_components(1) = 1;
            setappdata(tag_choice_panel,'n_components', n_components);     % store the number of components in this stimulus
                                                                  %  (array)
            current_component(1) = 1;
            setappdata(tag_choice_panel, 'current_component', current_component);
        else
            n_stimuli = size(entry.stimulus,2);
            setappdata(tag_choice_panel,'n_stimuli', n_stimuli);          % store the number of stimuli in this experiment
            if isappdata(tag_choice_panel,'current_stimulus')&&getappdata(tag_choice_panel,'current_stimulus')<=n_stimuli
                current_stimulus = getappdata(tag_choice_panel,'current_stimulus');
            else
                current_stimulus = 1;
                setappdata(tag_choice_panel,'current_stimulus', current_stimulus);
            end
            for idx_stimuli = 1:n_stimuli
                n_components(idx_stimuli) = size(entry.stimulus{idx_stimuli},2);
            end
            setappdata(tag_choice_panel,'n_components', n_components);     % store the number of components in this stimulus
                                                                           %  (array)
            if isappdata(tag_choice_panel,'current_component')
                current_component = getappdata(tag_choice_panel,'current_component');
            else
                current_component = size(n_components);
                setappdata(tag_choice_panel, 'current_component', current_component);
            end
        end
        
        % begin stimulus panel implementations
        if isappdata(tag_choice_panel, 'cur_panel_h')
            if getappdata(tag_choice_panel, 'cur_panel_h')~=get(hObject,'parent') % works if from the stimulus_number
                                                                                  % or stimulus_component_number panels
                stimulus_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.21, 0.01, 0.38,0.88],...
                    'BackgroundColor', [.13 .94, .55], 'title', '<stimulus> field options',...
                    'fontsize', 12, 'fontweight', 'bold');
                setappdata(tag_choice_panel, 'cur_panel_h', stimulus_panel);  % handle of this uipanel
                setappdata(tag_choice_panel, 'current_panel', 'stimulus_panel');
                setappdata(tag_choice_panel, 'cur_panel_saved', 0);
            else
                stimulus_panel = getappdata(tag_choice_panel, 'cur_panel_h');
            end
        else
            stimulus_panel = uipanel(mainfig, 'Units', 'normalized', 'Position', [0.21, 0.01, 0.38,0.88],...
                'BackgroundColor', [.13 .94, .55], 'title', '<stimulus> field options',...
                'fontsize', 12, 'fontweight', 'bold');
            setappdata(tag_choice_panel, 'cur_panel_h', stimulus_panel);  % handle of this uipanel
            setappdata(tag_choice_panel, 'current_panel', 'stimulus');
            setappdata(tag_choice_panel, 'cur_panel_saved', 0);
        end
        
    end
% --------------------/display_stimulus_fieldoptions----------------------

% display_experimental_animal_fieldoptions simply opens a panel with the experimental_animal
% field options
    function display_experimental_animal_fieldoptions(hObject, eventdata)
        a = 0;
    end
% --------------------/display_experimental_animal_fieldoptions----------------------

%% HELPER FUNCTIONS FOR SUBFIELDS

% screen_fieldoption_panel checks to ensure that you don't overwrite
% current info (1 = passed, 0 = failed)
    function result = screen_fieldoption_panel(hObject, eventdata)
        if isappdata(tag_choice_panel, 'current_panel')
            tag_choice_string = get(tag_choice_popup,'string');
            tag_choice_value = get(tag_choice_popup,'value');
            choice = tag_choice_string{tag_choice_value};
            if strmatch(choice, getappdata(tag_choice_panel, 'current_panel'))
                result = 1;
                return
            elseif isappdata(tag_choice_panel, 'cur_panel_saved')
                if getappdata(tag_choice_panel, 'cur_panel_saved')==0
                    choice = questdlg('You have not yet saved the data in this field.  Save now?',...
                        'Save Dialog', 'YES', 'NO', 'CLEAR DATA','NO');
                    switch choice
                        case 'YES'
                            save_entry(hObject, eventdata);
                            setappdata(tag_choice_panel, 'cur_panel_saved', 1)
                            result = 1;
                        case 'NO'
                            result = 0;
                        case 'CLEAR DATA'
                            result = 1;
                            delete(getappdata(tag_choice_panel, 'cur_panel_h'));
                    end
                else
                    result = 1;
                    delete(getappdata(tag_choice_panel, 'cur_panel_h'));
                end
            else
                result = 1;
            end
        else
        result = 1;
        end
    end
% -------------------/screen_fieldoption_panel----------------------------

% add_to_selected (takes the highlighted category from the available list handle and
% adds it to the selected_list).  returns the resultant selected list
    function add_to_selected(hObject, eventdata, available, selected, save_variable)
        values = get(available, 'value');
        strings = get(available, 'string');
        old_selected = get(selected, 'string');
        if ~iscell(old_selected)
           if old_selected == 'empty'
            old_selected = [];
           end
        end
        if ~isempty(old_selected)
            new_selected = [old_selected; strings(values)];
        else
            new_selected = strings(values);
        end
        set(selected, 'string', new_selected);
        selected = rot90(get(selected, 'string'));
        if size(strfind(save_variable, '.'),2)==1
            field = save_variable(1:strfind(save_variable,'.')-1);
            subfield = save_variable(strfind(save_variable,'.')+1:end);
            entry.(field).(subfield) = selected;
        else
            entry.(save_variable) = selected;
        end
        strings(values) = [];
        set(available, 'string', strings, 'value', 1);
    end
% -----------/add_to_selected----------------------------------------------

% remove_from_selected (takes the highlighted tag from the selected list and
% adds it to the available list.  retunrs the resultant selected list
    function remove_from_selected(hObject, eventdata, available, selected,save_variable)
        values = get(selected, 'value');
        strings = get(selected, 'string');
        old_selected = get(available, 'string');
        if ~iscell(strings)
           if strings == 'empty'
            return;
           end
        end
        if ~isempty(old_selected)
            new_selected = [old_selected; strings(values)];
        else
            new_selected = strings(values);
        end
        set(available, 'string', new_selected);
        strings(values) = [];
        set(selected, 'string', strings, 'value', 1);
        selected = rot90(get(selected, 'string'));
        if size(strfind(save_variable, '.'),2)==1
            field = save_variable(1:strfind(save_variable,'.')-1);
            subfield = save_variable(strfind(save_variable,'.')+1:end);
            entry.(field).(subfield) = selected;
        else
            entry.(save_variable) = selected;
        end
    end
% -----------/remove_from_selected----------------------------------------

% add_to_field function takes extra string arguments 'save_variable' which defines
% which subfield within the global <entry> value should be updated and
% 'type' which can be 'numeric' or 'string' (to allow for a check for
% validity before assigning the values. Varargin should be pairs of values.  Right now,
% this is designed to allow the user to identify a specific uicontrol handle to assign
% the new values into the 'string' category.  Designed with 'edit' boxes in mind
    function add_to_field(hObject, eventdata, save_variable, type, varargin)
        new_value = get(hObject, 'string');
        if isempty(new_value)
            return;
        end
        % checks to make sure entered data is of correct type
        switch type
            case 'numeric'
                if isempty(str2num(new_value))
                    errordlg('This field is numeric, please supply a numeric value')
                    set(hObject, 'string', '');
                    return;
                end
                new_value = str2num(new_value);
            case 'string'
                a = 0; % do nothing for now, just let the string be the string
        end
        % Update the existing field
        if size(strfind(save_variable, '.'),2)==1
            field = save_variable(1:strfind(save_variable,'.')-1);
            subfield = save_variable(strfind(save_variable,'.')+1:end);
            if isfield(entry,field)
                if isfield(entry.(field), subfield)
                    current_values = entry.(field).(subfield);
                else
                    current_values = [];
                end
            else
                current_values = [];
            end
            if isempty(current_values)
                entry.(field).(subfield) = new_value;
            else
                if ~ismember(entry.(field).(subfield), new_value)
                    entry.(field).(subfield)(end+1) = new_value;
                    if isnumeric(entry.(field).(subfield))
                        entry.(field).(subfield) = sort(entry.(field).(subfield));
                    end
                end
            end
            update_string = rot90(entry.(field).(subfield),3);
        else
            if isfield(entry, save_variable)
                current_values = entry.save_variable;
            else
                current_values = [];
            end
            if isempty(current_values)
                entry.(save_variable) = new_value;
            else
                entry.(save_variable)(end+1) = new_value;
            end
            update_string = rot90(entry.(save_variable),3);
        end
        
        % update the string
        if isnumeric(update_string)
            converted_string = cell(1);
            for idx = 1:size(update_string,1)
                converted_string(idx) = {num2str(update_string(idx))};
            end
            update_string = converted_string;
        end
        for idx_var = 1:size(varargin,2)
            if ishandle(varargin{idx_var})
                set(varargin{idx_var},'string',update_string,'max',size(update_string,1));
            end
        end
        set(hObject,'string', '');
    end
% --------------add_to_field-----------------------------------------------

% remove_from_field removes the selected item from the field
% (intended original use from popupmenu), requires variables
% 'save_variable' which is the subfield(s) under the <entry> structure
% and 'handle' which is the handle
% of the object that needs its current value removed
    function remove_from_field(hObject, eventdata, save_variable, handle)
        current_value = get(handle, 'value');
        current_string = get(handle, 'string');
        if size(strfind(save_variable, '.'),2)==1
            field = save_variable(1:strfind(save_variable,'.')-1);
            subfield = save_variable(strfind(save_variable,'.')+1:end);
            entry.(field).(subfield)(current_value) = [];
        else
            entry.(save_variable)(current_value) = [];
        end
        current_string(current_value) = [];
        set(handle,'string',current_string,'value',1);
    end
% ------/remove_from_field -----------------------------------------------

% change_value takes only (hObject,eventdata) and 'save_variable' which is
% a string that has the subfields (i.e. '.stimulus.conc_unit') within the
% <entry> structure.  Note: currently only works with string values
    function change_value(hObject, eventdata, save_variable)
       string = get(hObject, 'string');
       new_value = get(hObject, 'value');
       current = getappdata(get(hObject,'parent'),'current_idx');
       if size(strfind(save_variable, '.'),2)==1
            field = save_variable(1:strfind(save_variable,'.')-1);
            subfield = save_variable(strfind(save_variable,'.')+1:end);
            entry.(field){current}.(subfield) = string{new_value};
        else
            entry.(save_variable){current} = string{new_value};
        end
    end
% --------------/change_value---------------------------------------------

% update_stimulus
    function update_stimulus(hObject, eventdata)
        new_stimulus = get(hObject, 'value');
        if new_stimulus~=getappdata(get(hObject,'parent'), 'current_idx')
            setappdata(get(hObject,'parent'), 'current_idx', new_stimulus);
            display_stimulus_fieldoptions(hObject, eventdata);
        else
            return;
        end
    end
% -------------------/update_stimulus-------------------------------------

% add_stimulus
    function add_stimulus(hObject, eventdata)
        n_stimuli = size(entry.stimulus,2);
        entry.stimulus(n_stimuli+1) = {struct('category', {''}, 'identity', {''}, 'concentration', {[]}, ...
                'conc_unit', {''}, 'duration', {[]}, 'user_tag', {''}, 'strain', {''}, 'sex', {''})};
        setappdata(get(hObject, 'parent'),'current_idx',n_stimuli+1)
        display_stimulus_fieldoptions(hObject,eventdata);
    end
% --------------/add_stimulus---------------------------------------------

end