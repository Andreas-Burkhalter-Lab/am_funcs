function [entry] = cellcount_flat(xdb_id, varargin)
%CELLCOUNT_FLAT GUI interface for cell counting within regions of flattened images
%   This function takes in an argument in the form of an exp_db entry (.xdb
%   file) and opens up a GUI for interaction with the associated files
%   indicated within that .xdb entry
%
%   The inputs to the function are:
%       xdb_id: the unique identifier (same as base file name) of the experiment
%       varargin: only one current value is allowed, which is 'blind'.  If
%       user wishes to analyze expt blindly, he/she simply need enter
%       'blind' as a variable after the XDB_ID.
%
%   The inputs are used to open an experimental entry, upon which
%   additional fields will be defined by interaction with the GUI.
%
% DRAWING REGIONS:
%   By selecting an image and choosing the [DRAW REGION] button, user will
%       be able to define a subregion of interest within the image by clicking
%       on the image.  
%           DRAW CONTROLS: left-click: add point
%                          right-click: remove last point
%                          middle-click: finish region (complete the
%                                   circle) and disable drawing tool
%   After drawing the region, user will be prompted to enter a name 
%   for the subregion (i.e. 'glomerular, mitral').  The following fields
%   will be added to the saved .xdb structure upon choosing the [SAVE TO
%   XDB] button:
%       STRUCTURE FIELDS pertaining to REGIONS:
%           entry.hist_analysis.regions = coordinate points (x,y pairs)
%                                         defining the polygon drawn by the user
%                             regions = 1 x (number of images) cell with a 
%                                       cell array of length (number of
%                                       created regions).  Each individual
%                                       region is a cell array of
%                                       coordinate points
%           entry.hist_analysis.region_names = name of each defined region
%                             region_names = 1 x number of images cell,
%                                       each with a 1 x (number of regions
%                                       for this image) cell array of
%                                       string descriptors such as
%                                       'glomerular, mitral'
%           entry.hist_analysis.region_labels = the associated number of
%                                       each defined region (may be
%                                       eliminated in the future) same size
%                                       as region_names but with cell array
%                                       of integer labels
%         OPTIONAL FIELDS (added by clicking analysis buttons at bottom):
%           entry.hist_analysis.region_areas = the pixel area of each
%                                       user-defined region.  Same size as
%                                       region_names but with cell array of
%                                       doubles
%           entry.hist_analysis.um_per_pixels = # of microns per pixel for
%                                       images in this experiment.
%                                       Currently this is a 1-element cell
%                                       with a number 0 < x < 100.  Thus
%                                       all images are treated as having
%                                       the same number of um/pixel.
%
% SELECTING/COUNTING POSITIVE CELLS
%   Before or after saving regions, clicking the [COUNT CELLS] button
%           will enable user to click on the image to select cells (which will be marked by
%           a green circle.  The coordinates of each click will be stored in an
%           application variable for later use/saving.
%           entry.image_analysis.cell_positions = x, y coordinates of each cell identified
%       CELL COUNTING CONTROLS: 
%                          left-click: add point
%                          right-click: remove last point
%                          middle-click: finish counting cells and disable drawing tool
%
%       STRUCTURE FIELDS pertaining to CELL_POINTS:
%           entry.hist_analysis.cell_points = coordinate points (x,y pairs)
%                                         defining the cells identified the user
%         OPTIONAL FIELDS (set by choosing clicking analysis buttons)
%           entry.hist_analysis.cells_in_region = set by clicking [SORT
%                                         CELLS] button, counts the number
%                                         of cells which lie within each
%                                         defined region. 1 x (number of
%                                         images) cell array with a
%                                         1x(number of regions) number
%                                         array containing the cell counts
%
%   In case the user would like to use the information contained within the
%   saved .xdb information upon loading, this function returns the loaded
%   variables associated with the entry.  
%       PLEASE NOTE: this returned data does NOT change during
%       GUI use.  Only upon saving all input data and re-loading the GUI will
%       the user be able to view the structure contents of saved data.
%
%   Copyright 2007 Julian P. Meeks
%   Timothy Holy Lab
%
%   See also CALCULATE_REGION_AREA, SORT_CELLS_INTO_REGIONS, SAVE_XDB_ENTRY

%   Version History:
%       2007-08-28: Wrote outline and initial backbone
%       2007-09-06: Finished initial writing and documentation
%       2007-09-18: changed input params to include 'blind'
%                   eliminated need for username to be provided in function call

%% Load in entry, identify associated image files

% first, load in the entry
username = xdb_id(1:strfind(xdb_id, '_')-1);
entry = load(['/usr/lab/exp_db/' username filesep xdb_id '.xdb'], '-mat');
entry_copy = entry; % store original values for future reference and to reset if user decides

% second, determine if this person should be blind to the experimental
% date/filenames
if ~isempty(varargin)
    if varargin{1} == 'blind'
        blinded = 1;
    else
        blinded = 0;
    end
else 
    blinded = 0;
end

% next, find files which are saved in '.tif' format 
%   RESULT: is_tif is a matrix of 1s and 0s associated with file list in
%           entry.data_locations.  Basically a lookup table for which files
%           to attempt opening.
%       (TO-DO: add support for others, optimize searching?)
for idx_data_locs = 1:length(entry.data_locations)
    is_tif(idx_data_locs) = strwcmp('*.tif',entry.data_locations{idx_data_locs});
end
tif_indices = find(is_tif); % tif_indices(:) are the indices of entry.data_locations
                            % which have been verified as containing *.tif
                            % in filename (not necessarily verified images
                            % though)

current_tif = tif_indices(1);   % display first .tif file upon opening GUI

%% Setup GUI

% Initialize masterfig, hiding for now --------------------------------
masterfig = figure('Visible', 'on', 'Name', 'Cell Counting GUI', ...
                   'Position', [0 24 1200 700]); 
% --------------/masterfig---------------------------------------------

% Initialize main_image to contain first original image ---------------
image_axis = axes('Units', 'Pixels', 'Position', [210 5 985 690]);
set(image_axis, 'Units', 'Normalized');
% --------------/main_image--------------------------------------------

% Initialize image_hg (hggroup) to contain all items pertaining to the
% image so that they can be treated as a grouped object
 image_hg = hggroup('Parent', image_axis);
% --------------/image_hg----------------------------------------------

% Initialize xdb_id panel for top-left of UI ----------------------------
id_panel = uipanel(masterfig, 'Units', 'pixels', 'Position', [5 645 200 50]);
set(id_panel, 'Units', 'normalized', 'BackgroundColor', 'yellow');
id_label = uicontrol(id_panel, 'Units', 'normalized',...
               'Style', 'text', 'String', ['Experiment ID : ' entry.id], ...
               'Position', [0.01,0.61,0.98,0.38], 'FontSize', 12);
set(id_label, 'FontUnits', 'normalized', 'BackgroundColor', 'yellow');
if blinded == 0
    date_label = uicontrol(id_panel, 'Units', 'normalized',...
               'Style', 'text', 'String', ['Experiment date : ' datestr(datenum(entry.start_date))], ...
               'Position', [0.01,0.01,0.98,0.58], 'FontSize', 12);
    set(date_label, 'FontUnits', 'normalized', 'BackgroundColor', 'yellow');
end
% -------------/id_panel-----------------------------------------------

% ---------------------------------------------------------------------
% Initialize file_browse_panel to display which files are being displayed
%            and how many are defined in the current experiment
    % -- Global variables for this panel -----------
    %       already defined: tif_indices
    % define global variable tif_filename
    tif_filename = strip_path(entry.data_locations{current_tif}) ; % removes directory portion
    % ------------ /global variables ---------------
file_browse_panel = uipanel(masterfig, 'Units', 'pixels', 'Position', [5 580 200 60]);
set(file_browse_panel, 'Units', 'normalized', 'BackgroundColor', 'cyan');
% text box
file_browse_text_1 = uicontrol(file_browse_panel, 'Units', 'normalized',...
                               'Style', 'text', 'String', ' Image file ',...
                               'Position', [0.01,0.51,0.38,0.38], ...
                               'BackgroundColor', 'cyan', 'FontSize', 12);
set(file_browse_text_1, 'FontUnits', 'normalized');
% text box
file_browse_text_2 = uicontrol(file_browse_panel, 'Units', 'normalized',...
                               'Style', 'text', 'String', ['of ' num2str(length(tif_indices))],...
                               'Position', [0.55,0.51,0.30,0.38], ...
                               'BackgroundColor', 'cyan', 'FontSize', 12);
set(file_browse_text_1, 'FontUnits', 'normalized');
% popup menu
file_browse_popup = uicontrol(file_browse_panel, 'Units', 'normalized', ...
                                'Style', 'PopUpMenu', 'String', num2str(tif_indices(:)),...
                                'Position', [0.41,0.65,0.13,0.34],...
                                'FontSize', 12, 'Callback', @file_browse_update);
set(file_browse_popup, 'Value', 1, 'FontUnits', 'normalized');
% updating file name text box

if blinded == 0
    file_browse_filename = uicontrol(file_browse_panel, 'Units', 'normalized',...
                               'Style', 'text', 'String', tif_filename,...
                               'Position', [0.01,0.01,0.98,0.54], ...
                               'BackgroundColor', 'cyan', 'FontSize', 11);
    set(file_browse_text_1, 'FontUnits', 'normalized');
end
%----------/file_browse_panel-------------------------------------------
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% Initialize image_manip_panel to allow basic image manipulation for the
%   currently displayed figure
image_manip_panel = uipanel(masterfig, 'Units', 'pixels', 'Position', [5 475 200 100]);
set(image_manip_panel, 'Units', 'normalized', 'BackgroundColor', 'green');
% grayscale button
grayscale_btn = uicontrol(image_manip_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.01,0.66,0.28, 0.3], 'String', 'grayscale', ...
                          'Callback', @grayscale_image);
set(grayscale_btn, 'Units', 'normalized');
% jet button
jet_btn = uicontrol(image_manip_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.01,0.33,0.28, 0.3], 'String', 'jet', ...
                          'Callback', @jet_image);
% brightness slider
bright_slider = uicontrol(image_manip_panel, 'Style', 'slider', 'Units', 'normalized',...
                          'Position', [0.31,0.8,0.6,0.18], 'Value', 0, ...
                          'Max', 1, 'Min', -1, 'SliderStep', [0.01 0.1],...
                          'Callback', @bright_slide);
% brightness text
bright_text = uicontrol(image_manip_panel, 'Style', 'text', 'Units', 'normalized',...
                          'Position', [0.32,0.65,0.66,0.14], 'FontSize', 12,...
                          'String', ['Brightness : ' num2str(get(bright_slider, 'Value'))]);
set(bright_text, 'FontUnits', 'normalized');
% auto-contrast button
contrast_btn = uicontrol(image_manip_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.41,0.31,0.33, 0.29], 'String', 'Contrast', ...
                          'Callback', @contrast_image);
%--------------------/image_manip panel --------------------------------

%-----------------------------------------------------------------------
% Initialize draw_region panel
draw_region_panel = uipanel(masterfig, 'Units', 'pixels', 'Position', [5 370 200 100]);
set(draw_region_panel, 'Units', 'normalized', 'BackgroundColor', 'magenta');
% select_region_btn pushbutton
regions_in_play = zeros(size(tif_indices));
select_region_txt1 = uicontrol(draw_region_panel, 'style', 'text', 'Units', 'normalized',...
                          'Position', [0.01,0.85,0.48, 0.14], 'String', 'Region :',...
                          'BackgroundColor', 'magenta');
select_region_txt2 = uicontrol(draw_region_panel, 'style', 'text', 'Units', 'normalized',...
                          'Position', [0.31,0.55,0.12, 0.14], ...
                          'String', ['of ' num2str(regions_in_play(current_tif))],...
                          'BackgroundColor', 'magenta');
select_region_popup = uicontrol(draw_region_panel, 'Units', 'normalized', ...
                                'Style', 'PopUpMenu', 'String', ' ',...
                                'Position', [0.11,0.59,0.16,0.16],...
                                'Callback', @select_region_update);
set(file_browse_popup, 'Value', 1, 'FontUnits', 'normalized');
% draw_region_btn pushbutton
draw_region_btn = uicontrol(draw_region_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.01,0.01,0.48,0.48], 'String', 'Draw Region', ...
                          'Callback', @draw_region);
% draw_region_btn pushbutton
delete_region_btn = uicontrol(draw_region_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.51,0.01,0.48,0.48], 'String', 'Delete Region', ...
                          'Callback', @delete_region);
% draw_region_btn pushbutton
rename_region_btn = uicontrol(draw_region_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.51,0.51,0.48,0.48], 'String', 'Rename Region', ...
                          'Callback', @rename_region);
%------------------------/draw_region_panel-----------------------------

%-----------------------------------------------------------------------
% Initialize count_cells_panel
count_cells_panel = uipanel(masterfig, 'Units', 'pixels', 'Position', [5 270 200 99]);
set(count_cells_panel, 'Units', 'normalized', 'BackgroundColor', [.75 .5 .5]);
% draw_region_btn pushbutton
count_cells_btn = uicontrol(count_cells_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.01,0.51,0.48,0.48], 'String', 'Count Cells', ...
                          'Callback', @count_cells);
count_cells_txt1 = uicontrol(count_cells_panel, 'style', 'text', 'Units', 'normalized',...
                          'Position', [0.51,0.70,0.48, 0.16], 'String', '# cells in image :',...
                          'BackgroundColor', [.75 .5 .5]);
count_cells_number = uicontrol(count_cells_panel, 'style', 'text', 'Units', 'normalized',...
                          'Position', [0.51,0.55,0.48, 0.16], 'String', num2str(0),...
                          'BackgroundColor', [.75 .5 .5]);

%------------------------/draw_region_panel-----------------------------

%-----------------------------------------------------------------------
% Initialize file_io_panel  - where the saving happens
file_io_panel = uipanel(masterfig, 'Units', 'pixels', 'Position', [5 170 200 99]);
set(file_io_panel, 'Units', 'normalized', 'BackgroundColor', [.75 .75 .25]);
% draw_region_btn pushbutton
save_info_btn = uicontrol(file_io_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.1,0.1,0.8,0.8], 'String', 'Save to XDB', ...
                          'Callback', @save_info);
%------------------------/file_io_panel-----------------------------

%-----------------------------------------------------------------------
% Initialize function_panel  - where the common analysis functions are
% linked to current data for convenient analysis throughput
function_panel = uipanel(masterfig, 'Units', 'pixels', 'Position', [5 20 200 149]);
set(function_panel, 'Units', 'normalized', 'BackgroundColor', [1 .75 .75]);
% draw_region_btn pushbutton
calc_areas_btn = uicontrol(function_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.01,0.51,0.48,0.48], 'String', 'Calc Region Areas', ...
                          'Callback', @calculate_areas);
sort_cells_btn = uicontrol(function_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.01,0.01,0.48,0.48], 'String', 'Sort Cells', ...
                          'Callback', @sort_cells);
set_um_pixel_btn = uicontrol(function_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.51,0.51,0.48,0.48], 'String', 'Set um/pixel', ...
                          'Callback', @set_um_pixel);
set_um_pixel_btn = uicontrol(function_panel, 'style', 'pushbutton', 'Units', 'normalized',...
                          'Position', [0.51,0.01,0.48,0.48], 'String', 'ALL FUNCTIONS', ...
                          'Callback', @do_all_functions);%------------------------/draw_region_panel-----------------------------

%% Initialize important application data (appdata)

%  for all region variables, 'm' = maximum length of the cells for each region
%                            'n' = size(tif_indices)

% image_axis -> regions = cell(m, n).  Initially 1-by-n . contains the
%                       coordinate points associated with a drawn region
% image_axis -> n_regions = cell(n): integer # of FULL regions per tif
% image_axis -> region_names = cell(m,n).  Contains string names of regions
% image_axis -> region_lines = cell(m,n).  Contains handles of line objects
%                                           for each drawn region.
% image_axis -> region_labels = cell(m,n). Contains handles of text labels
%                                           for each drawn region
% image_axis -> cell_points = cell(n)   Contains all points for each image
% image_axis -> cell_point_handles = cell(n) one handle for each drawn point

setappdata(image_axis, 'regions', cell(size(tif_indices)));
setappdata(image_axis, 'n_regions', cell(size(tif_indices)));
setappdata(image_axis, 'region_names', cell(size(tif_indices)));
setappdata(image_axis, 'region_lines', cell(size(tif_indices)));
setappdata(image_axis, 'region_labels', cell(size(tif_indices)));
setappdata(image_axis, 'cell_points', cell(size(tif_indices)));
setappdata(image_axis, 'cell_point_handles', cell(size(tif_indices)));

% load any previously stored information from this experiment entry
if isfield(entry, 'hist_analysis')
    if isfield(entry.hist_analysis, 'regions');
        setappdata(image_axis, 'regions', entry.hist_analysis.regions);
        temp_n_regions = cell(size(entry.hist_analysis.regions));
        for idx = 1:size(entry.hist_analysis.regions,2)
            temp_n_regions(idx) = {size(entry.hist_analysis.regions{idx},2)};
        end
        setappdata(image_axis, 'n_regions', temp_n_regions);
        clear temp_n_regions;
    end
    if isfield(entry.hist_analysis, 'region_names')
        setappdata(image_axis, 'region_names', entry.hist_analysis.region_names);
    end
    if isfield(entry.hist_analysis, 'region_labels')
        setappdata(image_axis, 'region_labels', entry.hist_analysis.region_labels);
    end
    if isfield(entry.hist_analysis, 'cell_points')
        setappdata(image_axis, 'cell_points', entry.hist_analysis.cell_points);
    end
end

%% Display initial image 
    img_h = zeros(size(tif_indices));   % array to contain image handles   
    img = cell(size(tif_indices));          % .. to contain image data
    img{tif_indices(1)} = imread(entry.data_locations{tif_indices(1)}, 'tiff');
    img_h(tif_indices(1)) = imshow(img{tif_indices(1)}, 'Parent', image_axis);
    set(image_axis, 'clim',[1000 10000]);
    file_browse_update(file_browse_popup);
    set(image_axis, 'Visible', 'off');
    setappdata(image_axis, 'current_map', colormap);
% ----------------------------------------------------------------------

%% UI Function Callbacks

% function file_browse_update-------------------------------------------
%       : to update and display the selected image based on chosen 'Value' of file_browse_popup
    function file_browse_update(hObject, eventdata)
        current_choice = get(hObject, 'Value');
        if isappdata(image_axis, 'clim_vals')
            climvals = getappdata(image_axis, 'clim_vals');
            if max(climvals{current_tif}~=get(image_axis, 'Clim'))~=0
                climvals{current_tif} = get(image_axis, 'Clim');
                setappdata(image_axis, 'clim_vals', climvals);
            end
        else
            climvals = cell(size(tif_indices));
            initial_clims = [1000 10000];
            for i = 1:size(tif_indices, 2)
                climvals{i} = initial_clims;
            end
            setappdata(image_axis,'clim_vals', climvals);
        end
        apply_clim = climvals{current_choice};
   hold off;
        img{tif_indices(current_choice)} = imread(entry.data_locations{tif_indices(current_choice)}, 'tiff');
        img_h(tif_indices(current_choice)) = imshow(img{tif_indices(current_choice)},...
                                                 'Parent', image_axis);
        set(image_axis, 'clim',[1000 10000]);
        if ~isempty(apply_clim)
            set(image_axis, 'Clim', apply_clim);
        end
        if strmatch('jet', getappdata(image_axis, 'map_flag'))
            colormap(image_axis,jet);
        elseif strmatch('gray',getappdata(image_axis, 'map_flag'))
            colormap(image_axis,gray);
        end
        set(image_axis, 'Visible', 'off');
        current_tif = current_choice;
        tif_filename = strip_path(entry.data_locations{current_tif});
        if blinded == 0;
            set(file_browse_filename, 'String', tif_filename);
        end
        
        % draw existing region lines
        regions = getappdata(image_axis, 'regions');
        region_lines = getappdata(image_axis, 'region_lines');
        region_names = getappdata(image_axis, 'region_names');
        n_regions = getappdata(image_axis, 'n_regions');
        for idx = 1:size(regions{current_tif},2)
           if ~isempty(regions{current_tif}{idx})
               draw_region_on_screen(idx, current_tif);
           end
        end
        if size(regions{current_tif},2)~=0
            % set the popup menu to list the correct number of areas
            set(select_region_popup, 'String', num2str(rot90(1:n_regions{current_tif},3)));
            set(select_region_popup, 'Value', n_regions{current_tif});
            set(select_region_txt2, 'String', ['of ' num2str(n_regions{current_tif})]);
        else
            set(select_region_popup, 'String', ' ');
            set(select_region_popup, 'Value', 1);
            set(select_region_txt2, 'String', ['of ' num2str(0)]);
        end
        
        % draw any existing cell_points
        draw_cellpoints_on_screen(current_tif, 'redraw');
        update_cellpoint_number(current_tif);
     hold on;
    end
%----------/file_browse_update------------------------------------------

% function grayscale_image (does as it says)
    function grayscale_image(hObject, eventdata)
        colormap(image_axis, gray);
        if ~isappdata(image_axis, 'map_flag')
            setappdata(image_axis, 'map_flag', 'gray');
        end
        if ~isappdata(image_axis, 'orig_graymap')
            setappdata(image_axis, 'orig_graymap', colormap(image_axis));
        end
        setappdata(image_axis, 'map_flag', 'gray');
        set(bright_slider, 'Value', 0);
        bright_slide(bright_slider);
    end
%----------/grayscale_image---------------------------------------------

% function jet_image (does as it says)
    function jet_image(hObject, eventdata)
        colormap(image_axis, jet);
        if ~isappdata(image_axis, 'map_flag')
            setappdata(image_axis, 'map_flag', 'jet');
        end
        if ~isappdata(image_axis, 'orig_jetmap')
            setappdata(image_axis, 'orig_jetmap', colormap(image_axis));
        end
        setappdata(image_axis, 'map_flag', 'jet');
        set(bright_slider, 'Value', 0);
        bright_slide(bright_slider);
    end
%----------/jet_image---------------------------------------------------

% function bright_slide brightens/darkens the image
    function bright_slide(hObject,eventdata)
        bright_val = get(hObject, 'Value');
        map_flag = getappdata(image_axis, 'map_flag');
        if strmatch('gray', map_flag)
            orig_map = getappdata(image_axis, 'orig_graymap');
        elseif strmatch('jet', map_flag)
            orig_map = getappdata(image_axis, 'orig_jetmap');
        end
        cur_map = orig_map+bright_val;
        cur_map(cur_map>1)=1;
        cur_map(cur_map<0)=0;
        setappdata(image_axis, 'current_map', cur_map);
        colormap(image_axis, getappdata(image_axis, 'current_map'));
        set(bright_text, 'String', ['Brightness : ' num2str(bright_val)]);
    end            
%-------------/bright_slide---------------------------------------------
        
% function contrast_image uses contrast_btn to perform autocontrast
    function contrast_image(hObject, eventdata)
        imcontrast(image_axis);
    end
%-----------------/contrast_image---------------------------------------

% function draw_region : to allow the user to select a -----------------
    function draw_region(hObject, eventdata)
        if ~isappdata(image_axis, 'regions')
            error('The application variable ''regions'' does not exist.');
        end

    % wipe any existing button functions in the current figure:
        set(get_parent_fig(image_axis), 'WindowButtonUpFcn', []);
        set(get_parent_fig(image_axis), 'windowButtonDownFcn', []);
    
    % set button 'down' function to call make_region_point( )
        install_mouse_event_handler(image_axis, 'down', @make_region_point);
    end
%---------------/draw_region--------------------------------------------


% function make_region_point to log the position of each button click on
%                       the image of interest (displayed in image_axis)
% 
    function result = make_region_point(sender, eventdata)
        
        regions = getappdata(image_axis, 'regions'); % retrieve current region data
        n_regions = getappdata(image_axis, 'n_regions');
        region_names = getappdata(image_axis, 'region_names');
        
        % assert: under normal operation, these regions will exist for
        % this function to be called (no specific error check right now though)
        
        if isempty(n_regions{current_tif})
            cur_region = 1;
        else
            cur_region = n_regions{current_tif}+1;
        end
        
        if cur_region>size(regions{current_tif},2)
            new_region = [];
        else
            new_region = regions{current_tif}{cur_region};  % local copy of current_tif's region(s)
        end
        if is_button_down(sender, 'left')
            pixel = get(sender, 'currentpoint');
            if isempty(new_region)
                new_region = {pixel(1,1:2)};
            else
                new_region(end+1) = {pixel(1, 1:2)};
            end
        elseif is_button_down(sender, 'right')
            if ~isempty(new_region)
                last = size(new_region,2);
                while isempty(new_region{last})
                    last = last-1;
                end
                   new_region(last) = [];
            end
        elseif is_button_down(sender, 'middle')
            if ~isempty(new_region)
                if size(new_region,2)>2
                    if isempty(new_region(end))
                        new_region(end) = new_region(1);
                    else
                        new_region(end+1) = new_region(1);
                    end
                    disp('ended region successfully');
                    n_regions(current_tif) = {cur_region};
                else
                    new_region = [];
                    errormsg = 'You must draw regions with at least 3 points!';
                end
                
                % set the current region name via user input:
                region_names{current_tif}{cur_region} = inputdlg('Enter the region name');
                setappdata(image_axis, 'region_names', region_names);
                % disable the region maker until another button press
                uninstall_mouse_event_handler(image_axis, 'down', @make_region_point);
                set(get_parent_fig(image_axis), 'WindowButtonUpFcn', []);
                set(get_parent_fig(image_axis), 'windowButtonDownFcn', []);
            end
            if size(regions{current_tif},2)~=0
                % set the popup menu to list the correct number of areas
                set(select_region_popup, 'String', num2str(rot90(1:n_regions{current_tif},3)));
                set(select_region_popup, 'Value', n_regions{current_tif});
                set(select_region_txt2, 'String', ['of ' num2str(n_regions{current_tif})]);
            else
               set(select_region_popup, 'String', ' ');
               set(select_region_popup, 'Value', 1);
               set(select_region_txt2, 'String', ['of ' num2str(0)]);
            end
        end
        regions{current_tif}{cur_region} = new_region;
        setappdata(image_axis, 'regions', regions);
        setappdata(image_axis, 'n_regions', n_regions);

        % assert: click has caused a new point to be added, last point to
        % be deleted, or termination of a region.
   
        draw_region_on_screen(cur_region, current_tif)
        highlight_region(cur_region, current_tif);
        if exist('errormsg', 'var')
            result = errormsg;
        else
            result = 0;
        end
    end
% -----------------/make_region_point ----------------------------------

% function draw_region_on_screen to display the current points in a given
% region on the current tif image.
    function handle = draw_region_on_screen(cur_region, current_tif)
        region_lines = getappdata(image_axis, 'region_lines');
        regions = getappdata(image_axis, 'regions');
        region_names = getappdata(image_axis, 'region_names');
        region_labels = getappdata(image_axis, 'region_labels');
        cur_points = regions{current_tif}{cur_region};
       
        % separate x, y coords for line drawing
        if ~isempty(cur_points)
            for idx = 1:length(cur_points)
                x_pts(idx) = cur_points{idx}(1);
                y_pts(idx) = cur_points{idx}(2);
            end

            if ~isempty(region_names{current_tif})&&size(region_names{current_tif},2)>=cur_region
                min_x = min(x_pts); max_y = max(y_pts);     % store values for text marker
                text_label = strcat('[', num2str(cur_region), ']:', region_names{current_tif}{cur_region});
            end
            
            if (size(region_lines{current_tif},2)>=cur_region)&&ishandle(region_lines{current_tif}{cur_region})
                delete(region_lines{current_tif}{cur_region});
            end
            region_lines{current_tif}{cur_region} = line(x_pts,y_pts,...
                                  'Parent',image_axis,'Color', 'red', 'LineWidth', 2,...
                                  'Marker', 'none');
            if exist('text_label', 'var')
                region_labels{current_tif}{cur_region} = text(min_x, max_y+20, text_label,...
                                                        'Color', 'red', 'FontSize', 12);
            end
        end
        
        setappdata(image_axis, 'region_lines', region_lines);
        setappdata(image_axis, 'region_labels', region_labels);
    end

% function select_region_update to allow user to select one of the regions
% currently assigned for the current tif: popupmenu callback
    function result = select_region_update(sender, eventdata)
        cur_region = get(sender, 'Value');
        highlight_region(cur_region, current_tif);
        result = 0;
    end
% ----------/select_region_update---------------------------------------

% function highlight_region to take the selected region and turn it
% white/bold
    function result = highlight_region(cur_region, current_tif)
        region_lines = getappdata(image_axis, 'region_lines');
        region_labels = getappdata(image_axis, 'region_labels');
        for idx = 1:size(region_lines{current_tif},2)
            if ~isempty(region_lines{current_tif}{idx})
                set(region_lines{current_tif}{idx}, 'Color', 'red', 'LineWidth', 2);
            end
        end
        for idx = 1:size(region_labels{current_tif},2)
            if ~isempty(region_labels{current_tif}{idx})
                set(region_labels{current_tif}{idx}, 'Color', 'red', 'FontSize', 12);
            end
        end
        set(region_lines{current_tif}{cur_region}, 'Color', 'white', 'LineWidth', 2);

        if ~isempty(region_labels{current_tif})&&size(region_labels{current_tif},2)>=cur_region
            set(region_labels{current_tif}{cur_region}, 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
        end
        
        result = 0;
    end

% function delete_region : to allow the user to select a ---------------
    function delete_region(hObject, eventdata)
      hold off 
      regions = getappdata(image_axis, 'regions');
      region_lines = getappdata(image_axis, 'region_lines');
      region_names = getappdata(image_axis, 'region_names');
      n_regions = getappdata(image_axis, 'n_regions');
      region_labels = getappdata(image_axis, 'region_labels');
           
      selected_region = get(select_region_popup, 'Value');
           
      %%% put in a check for the user to verify
           
      regions{current_tif}(selected_region) = [];
      
      set(region_lines{current_tif}{selected_region}, 'Visible', 'off');
      delete(region_lines{current_tif}{selected_region});
      region_lines{current_tif}(selected_region) = [];
      set(region_labels{current_tif}{selected_region}, 'Visible', 'off');
      delete(region_labels{current_tif}{selected_region});
      region_labels{current_tif}(selected_region) = [];
      
      region_names{current_tif}(selected_region) = [];
      if n_regions{current_tif}>0
          n_regions(current_tif) = {n_regions{current_tif}-1};
      else
          n_regions(current_tif) = {0};
      end
      
      setappdata(image_axis, 'regions', regions)
      setappdata(image_axis, 'region_names', region_names);
      setappdata(image_axis, 'region_lines', region_lines);
      setappdata(image_axis, 'n_regions', n_regions);
      setappdata(image_axis, 'region_labels', region_labels);
           
      for idx = 1:size(regions{current_tif},2)
          draw_region_on_screen(idx, current_tif)
      end

      if size(regions{current_tif},2)~=0
          % set the popup menu to list the correct number of areas
          set(select_region_popup, 'String', num2str(rot90(1:n_regions{current_tif},3)));
          set(select_region_popup, 'Value', n_regions{current_tif});
          set(select_region_txt2, 'String', ['of ' num2str(n_regions{current_tif})]);
      else
          set(select_region_popup, 'String', ' ');
          set(select_region_popup, 'Value', 1);
          set(select_region_txt2, 'String', ['of ' num2str(0)]);
      end
      hold on;
    end
%---------------/delete_region------------------------------------------

% function rename_region : to allow the user to rename a region 
    function result = rename_region(hObject, eventdata)
        region_names = getappdata(image_axis, 'region_names');
        selected_region = get(select_region_popup, 'Value');
        new_region_name = inputdlg('Rename region: ', 'Rename region dialog',...
                                   1, region_names{current_tif}{selected_region});
        region_names{current_tif}{selected_region} = new_region_name;
        result = 0;
        setappdata(image_axis, 'region_names', region_names);
        draw_region_on_screen(selected_region, current_tif);
    end
%---------------/rename_region------------------------------------------

% function count_cells : to allow the user to pick cells using the mouse 
    function result = count_cells(hObject, eventdata)
    % wipe any existing button functions in the current figure:
        set(get_parent_fig(image_axis), 'WindowButtonUpFcn', []);
        set(get_parent_fig(image_axis), 'windowButtonDownFcn', []);
    
    % set button 'down' function to call make_region_point( )
        install_mouse_event_handler(image_axis, 'down', @make_cell_point);
    end
%---------------/count_cells--------------------------------------------

% function make_cell_point: mouse picker: ------------------------------
%   left button: new point
%   right button: delete last point
%   middle point: finish counting
    function result = make_cell_point(sender, eventdata)
        cell_points = getappdata(image_axis, 'cell_points');
        cell_point_handles = getappdata(image_axis, 'cell_point_handles');
        if isempty(cell_points{current_tif})
            current_points = [];
        else
            current_points = cell_points{current_tif};
        end
        
        % button choice portion --------------------------------
        if is_button_down(sender, 'left')
            pixel = get(sender, 'currentpoint');
            if isempty(current_points)
                current_points = {pixel(1,1:2)};
            else
                current_points(end+1) = {pixel(1, 1:2)};
            end
        elseif is_button_down(sender, 'right')
            if ~isempty(current_points)
                last = size(current_points,2);
                while isempty(current_points{last})
                    last = last-1;
                end
                   current_points(last) = [];
                   if size(cell_point_handles{current_tif},2)>=last
                       delete(cell_point_handles{current_tif}{last})
                   end
            end
        elseif is_button_down(sender, 'middle')
            if ~isempty(current_points)
                disp('done counting points');
                setappdata(image_axis, 'cell_points', cell_points);
            end
            % wipe any existing button functions in the current figure:
            uninstall_mouse_event_handler(image_axis, 'down', @make_cell_point);
            set(get_parent_fig(image_axis), 'WindowButtonUpFcn', []);
            set(get_parent_fig(image_axis), 'windowButtonDownFcn', []);    
        end % end button choice --------------------------------
        
        % save data to the application
        cell_points(current_tif) = {current_points};
        setappdata(image_axis, 'cell_points', cell_points);
        
        % plot the points
        draw_cellpoints_on_screen(current_tif, 'new');
        update_cellpoint_number(current_tif);
        
        result = 0;
    end
%---------------/make_cell_point----------------------------------------

% function draw_cellpoints_on_screen to draw chosen cells for the current tif
    function draw_cellpoints_on_screen(current_tif, toggle)
        % toggle values: 'new', 'redraw'
        cell_points = getappdata(image_axis, 'cell_points');
        cell_point_handles = getappdata(image_axis, 'cell_point_handles');
        
        set(image_axis, 'NextPlot', 'add'); % ensure "PLOT" command doesn't overwrite image
    switch toggle
        case 'redraw'
            for idx = 1:size(cell_points{current_tif}, 2);
                x_pts(idx) = cell_points{current_tif}{idx}(1);
                y_pts(idx) = cell_points{current_tif}{idx}(2);
                cell_point_handles{current_tif}{idx} = ...
                        plot(x_pts(idx), y_pts(idx),'go',...
                        'Parent', image_axis);
            end
        case 'new'
            hold on;
            plot(cell_points{current_tif}{end}(1), cell_points{current_tif}{end}(2),...
                    'go', 'Parent', image_axis);
            hold off;
    end
        setappdata(image_axis, 'cell_point_handles', cell_point_handles);    
    end
%---------------/draw_cellpoints_on_screen------------------------------

% function update_cellpoint_number: to update the text number of cells
%                                   picked by the user in the current tif
    function update_cellpoint_number(current_tif)
        cell_points = getappdata(image_axis, 'cell_points');
        n_cell_points = size(cell_points{current_tif},2);
        set(count_cells_number, 'String', num2str(n_cell_points));
    end
%-------------/update_cellpoint_number----------------------------------

% function save_info: to save new fields to 'entry'
%
    function result = save_info(hObject, eventdata)
        if ~isfield(entry, 'hist_analysis')
            entry.analysis = struct;
        end
        entry.hist_analysis.regions = getappdata(image_axis, 'regions');
        entry.hist_analysis.region_names = getappdata(image_axis, 'region_names');
        entry.hist_analysis.cell_points = getappdata(image_axis, 'cell_points');
        entry.hist_analysis.region_labels = getappdata(image_axis, 'region_labels');
        
        % if region_areas calculated, ensure they are in the correct size before
        % saving
        if isappdata(image_axis, 'region_areas')
            verified = 0;
            region_areas = getappdata(image_axis, 'region_areas');
            if size(region_areas,2) == size(entry.hist_analysis.regions, 2)
                for idx_img = 1:size(region_areas,2)
                    n_region_areas(idx_img) = size(region_areas{idx_img},2);
                    n_regions(idx_img) = size(entry.hist_analysis.regions{idx_img},2);
                end
                if n_region_areas == n_regions
                    verified = 1;
                end
            end
            if verified == 1
                entry.hist_analysis.region_areas = region_areas;
            end
        end  % end region_areas validation
       
        % if cells_in_regions calculated, ensure they are in the correct size before
        % saving
        if isappdata(image_axis, 'cells_in_regions')
            verified = 0;
            cells_in_regions = getappdata(image_axis, 'cells_in_regions');
            if size(cells_in_regions,2) == size(entry.hist_analysis.regions, 2)
                for idx_img = 1:size(cells_in_regions,2)
                    n_c_i_r(idx_img) = size(cells_in_regions{idx_img},2);
                    n_regions(idx_img) = size(entry.hist_analysis.regions{idx_img},2);
                end
                if n_c_i_r == n_regions
                    verified = 1;
                end
            end
            if verified == 1
                entry.hist_analysis.cells_in_regions = cells_in_regions;
            end
        end  % end region_areas validation        
        
        save_xdb_entry(entry, struct('id', entry.id));
    end
%--------------------/save_info-----------------------------------------

% function calculate_areas: linked to calc_areas_btn, takes regions and
% passes them to calculate_region_area(regions)
    function calculate_areas(hObject, eventdata)
        regions = getappdata(image_axis, 'regions');
        region_areas = calculate_region_area(regions);
        setappdata(image_axis, 'region_areas', region_areas);
        save_info(hObject, eventdata);  % should only save if region_areas
                                        % are correct size
    end
%------------------/calculate_areas-------------------------------------

% function sort_cells: linked to sort_cells_btn, takes cell_points and
% regions and passes them to sort_cells_into_regions function
    function sort_cells(hObject, eventdata)
        regions = getappdata(image_axis, 'regions');
        cell_points = getappdata(image_axis, 'cell_points');
        cells_in_regions = sort_cells_into_regions(cell_points, regions);
        setappdata(image_axis, 'cells_in_regions', cells_in_regions);
        save_info(hObject, eventdata);
    end
%-------------------/sort_cells-----------------------------------------

% function set_um_pixel: allows user to input the number of microns/pixel
% (x,y dimensions assumed equal)
    function set_um_pixel(hObject, eventdata)
        um_per_pixel = inputdlg('What is the um/pixel conversion? :');
        if isstr(um_per_pixel{1}); um_per_pixel{1}=str2num(um_per_pixel{1});end;
        while (um_per_pixel{1} < 0 || um_per_pixel{1} > 100)
            um_per_pixel = inputdlg('Sorry, range must be 0 < x < 100 : ');
        end
        entry.hist_analysis.um_per_pixel = um_per_pixel;
        save_info(hObject,eventdata);
    end
%---------------------/set_um_pixel-------------------------------------

% function do_all_functions : choose all 3 analysis functions 
    function do_all_functions(hObject, eventdata)
        calculate_areas(hObject, eventdata);
        sort_cells(hObject, eventdata);
        set_um_pixel(hObject, eventdata);
        save_info(hObject, eventdata);
    end
%-----------------/do_all_functions-------------------------------------

end