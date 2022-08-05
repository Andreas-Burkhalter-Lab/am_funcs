% Template for adding an entry to the exerimental database ('xdb')
% Just follow step-by-step!

%% 1) Standard set up that you don't have to change

% initialize a structure that you will add information to as you go
entry = struct;
entry.version = 1; 
    % this just means this is the first version of the database

%% 2) There are some basic fields that everyone is required to include:

% Who are you?
entry.investigator = 'julian';
    % If you aren't already on the list of investigators, you may have to
    % open the function "validate_xdb_entry_var.m" and add yourself to the
    % list in the first section, labeled "**** What investigators are
    % there?".
    %
    % You can see who's on the list already (and how their names are
    % spelled...) by typing "display_xdb_entry_options('investigator')" at
    % the command line and hitting return.

% What day did you start the experiment on?
year = [2007];  % put the year as a 4 digit number
month = [08]; % put the month as a 2 digit number
day = [23]; % put the day as a 2 digit number
entry.start_date = [year month day];

% What kind of experiment did you do?
tags = {'histology','fluoro_jade','vno_aob_hemihead','section', 'confocal_imaging'};
    % Fill this field in with things like 'single_electrode' or 'singing'.
    % You can put as many tags as are relevant to your experiment - the
    % more the better!
    %
    % Enter them as strings (test with single quote marks around it) - 
    % for example, you might end up with something like:
    %     tags = {'single_electride','aob','sorted'};
    %
    % For a complete list of tags that are currently allowed, type
    % "display_xdb_entry_options('tags')" at the command line
    % and hit return.
entry.tags = tags;
    
% Where is your data for this experiment stored?
entry.data_locations = {strcat(pwd,filesep,'2007_08_23_1hrbench_section2_avg_merge.tif'),...
                       strcat(pwd,filesep,'2007_08_23_1hrbench_section3_avg_merge.tif'),...
                       strcat(pwd,filesep,'2007_08_23_1hrbench_section4_avg_merge.tif'),...
                       strcat(pwd,filesep,'2007_08_23_1hrbench_section5_avg_merge.tif'),...
                       };
    
    % example: entry.datalocations =
    % {'/home/rebecca/matlabdata/experiemntaldata/2007_08_21/experiment1/'}
    %
    % If your data is in more than one place, you can include both of them:
    % entry.datalocations = {'/home/rebecca/foldername/date/exp1/',...
    % 'home/rebecca/foldername/date/exp2/'};
    
%% Check for existing experiment using cookie file
if ~isempty(dirbyname('*.cookie'))
    existing_entry = search_xdb(@query_callback_check_datafiles);
else
    existing_entry = [];
end
if ~isempty(existing_entry)
    id_list = cellfun(@get_xdb_ids_from_entry, existing_entry, 'UniformOutput', false);
    prompt_string = strcat('Exp_db entry or entries \n', sprintf('%s\n',id_list{:}), '\nexist. Overwrite .xdb file? (y/n) : ');
    if input(prompt_string, 's')=='y'
        if length(existing_entry)>1
            overwrite_old_xdb = 'true';
            prompt_string = 'Which file in the above list should be overwritten? (numeric top to bottom) : ';
            file_choice = input(prompt_string);
        else
            overwrite_old_xdb = 'true';
            file_choice = 1;
        end
        for idx_data_locs = 1:length(entry.data_locations(:))
            if isempty(strmatch(entry.data_locations{idx_data_locs}, existing_entry{file_choice}.data_locations(:)))
                existing_entry{file_choice}.data_locations = cat(2, existing_entry{file_choice}.data_locations, entry.data_locations{idx_data_locs});
            end
        end
        entry = existing_entry{file_choice};
    else
        overwrite_old_xdb = 'false';
    end
end
    
    
%% Optional Fields

% EXPERIMENTAL_ANIMAL (pertaining to animal being experimented upon)
if isfield(entry, 'experimental_animal')        % if field exists, prompt to overwrite
    prompt_string = '"Experimental Animal" field exists, update/replace values? (y/n) : ';
    exp_animal_choice = input(prompt_string, 's');
end
if ~isfield(entry,'experimental_animal')||(exp_animal_choice=='y')  % if new field or chose to overwrite
        entry.experimental_animal.age = 69;     % How many days old is the animal?
        entry.experimental_animal.sex = 'M';    % what sex is the animal? ('M' or 'F')
        entry.experimental_animal.strain = 'B6D2F1/J';  % which strain is the animal?
end


% SECTION (pertaining to tissue sectioning procedure)
entry.section.units = 'um';                     % section thickness units
entry.section.thickness = 75;                   % section thickness
entry.section.orientation = 'sagittal';         % i.e. 'sagittal', 'coronal', 'horizontal', etc.
entry.section.instrument = 'Vibratome 1000';    % which sectioning equipment?

% HISTOLOGY (pertaining to histological processing)
entry.histology.method = 'fluoro_jade_c';   % i.e. 'fluoro_jade_c', 'nissl', 'hematoxylin_and_eosin', etc.
entry.histology.mounted = 'true';           % Was tissue on slides for staining (false if free-floating)
entry.histology.success = 'true';           % if 'false' analyze with extreme caution

% VNO_AOB_HEMIHEAD (pertaining to the hemi-head preparation)
entry.vno_aob_hemihead.score = 9;       % subjective score (1 = awful, 10 = perfect) (<5 analyze with caution)
entry.vno_aob_hemihead.chamber_version = 2.0;   % which model of perfusion chamber used?
entry.vno_aob_hemihead.flowrate = 7.5;   % how fast was flow applied? (in mL/min)
entry.vno_aob_hemihead.temperature = 33;    % what was bath temperature during experiment? (degrees C)
entry.vno_aob_hemihead.cannula_size = NaN; % if applicable, what diameter cannula? (in um)

% CONFOCAL_IMAGING (pertaining to settings used for confocal microscopy)
entry.confocal_imaging.microscope = 'FV500b';   % Which microscope used?
entry.confocal_imaging.dye = 'FITC';    % Which pre-set dye used (if applicable)
entry.confocal_imaging.laser_percent = 1;  % what % laser applied?
entry.confocal_imaging.z_step = 2.5;    % how many microns between z_steps?
entry.confocal_imaging.format = 'tiff'; % which image format?
entry.confocal_imaging.bits = 8;        % how many bits per pixel?
entry.confocal_imaging.pixel_dimensions = [1024 1024];  % original image size (format [x y])?
entry.confocal_imaging.flatten_method = 'average';  % if image flattened, by what method (from MetaMorph typically)
entry.confocal_imaging.montaged = 'true';   % if composite (montage) image created, set to 'true'

% FLUORO_JADE (pertaining to fluoro_jade experiments)
entry.fluoro_jade.expt_timepoint = 60;          % time in minutes
entry.fluoro_jade.expt_condition = 'bench';

%% Validate and save!

% put walk through of validation, and goign through the error_messages that
% you get, and how to fix common problems...
errors = validate_xdb_entry_var(entry);
if errors ~= 0
    errors
    break;
end

% save step: picking a name and saving
% Evaluate user decisions for dealing with existing .xdb entries and choose
% appropriate save step
if (isfield(entry,'id'))&&exist('overwrite_old_xdb', 'var')
    switch(overwrite_old_xdb)
        case 'true'
        options.id = entry.id;                              % if user has been prompted and chose to overwrite,
        saved_entry_id = save_xdb_entry(entry, options);    % write new .xdb file using existing ID #
        case 'false'
            if (exist(entry.associated_ids))                % if user chooses not to overwrite, proceed
            entry.associated_ids = cat(2,associated_ids,entry.id);  % by choosing a new experiment.id
            else                                            % but noting the overlap
            entry.associated_ids = {entry.id};              % in field entry.associated_ids
            end
            saved_entry_id = save_xdb_entry(entry);
    end
elseif isfield(entry,'id')
    prompt_string = sprintf('An existing entry id (%s) was found, overwrite? (y/n) : ', entry.id);
    if input(prompt_string, 's')=='y'                   % if user was NOT prompted previously
        options.id = entry.id;                          % but existing id was assigned, prompt now
        saved_entry_id = save_xdb_entry(entry, options);
    else
        if (exist(entry.associated_ids))                % if user chooses not to overwrite, proceed
            entry.associated_ids = cat(2,associated_ids,entry.id);  % by choosing a new experiment.id
        else                                            % but noting the overlap
            entry.associated_ids = {entry.id};          % in field entry.associated_ids
        end
        saved_entry_id = save_xdb_entry(entry);         
    end
elseif exist('existing_entry', 'var')&&~isempty(existing_entry)
    % assert: user did not indicate overwrite but existing entry(ies) exist,
    % Therefore create a new .xdb file with a new unique id # and links to
    % existing in entry.associated_ids
    if exist('existing_entry','var')&&isfield(existing_entry{:},'associated_ids')     % if user chooses not to overwrite, proceed
            entry.associated_ids = cat(2,associated_ids,existing_entry{:}.id);  % by choosing a new experiment.id
    elseif exist('existing_entry','var')&&isfield(existing_entry{:},'id')            % but noting the overlap
            entry.associated_ids = {existing_entry{:}.id};                    % in field entry.associated_ids
    end        
    saved_entry_id = save_xdb_entry(entry);
else
    % assert: user did not indicate overwrite, and no existing or new
    % entry.ids exist, therefore proceed to write a new unique entry with
    % no associated_ids
    saved_entry_id = save_xdb_entry(entry);
end
    

% save local cookie containing strings for 'entry.id' and the data
% locations
cookie_filename = '2007_08_23_bench.cookie';
if exist(cookie_filename,'file')&&exist('overwrite_old_xdb','var')
    switch overwrite_old_xdb
        case 'true'
            fprintf(1, '\nNote: existing cookie file may be updated along with .xdb file!\n\n');
        case 'false'
            while exist(cookie_filename, 'file')
               prompt_string = sprintf('Cookie file %s already exists. Please indicate new filename (*.cookie) : ', cookie_filename);
               cookie_filename = input(prompt_string,'s');
            end
    end
end
cookie.filename = strcat(pwd,filesep,cookie_filename);
cookie.xdb_id = strcat('/usr/lab/exp_db/',entry.investigator,filesep, saved_entry_id,'.xdb');
cookie.data_locations = entry.data_locations;

save_xdb_cookie(cookie);


% mention where to go when you want to acutally search for something...