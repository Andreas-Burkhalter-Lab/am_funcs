% Template for adding an entry to the exerimental database ('xdb')
% Just follow step-by-step!

%% 1) Standard set up that you don't have to change

% initialize a structure that you will add information to as you go
entry = struct;
entry.version = 1; 
    % this just means this is the first version of the database

%% 2) There are some basic fields that everyone is required to include:

% Who are you?
entry.investigator = 'put_your_name_here';
    % If you aren't already on the list of investigators, you may have to
    % open the function "validate_xdb_entry_var.m" and add yourself to the
    % list in the first section, labeled "**** What investigators are
    % there?".
    %
    % You can see who's on the list already (and how their names are
    % spelled...) by typing "display_xdb_entry_options('investigator')" at
    % the command line and hitting return.

% What day did you start the experiment on?
year = [];  % put the year as a 4 digit number
month = []; % put the month as a 2 digit number
day = []; % put the day as a 2 digit number
entry.start_date = [year month day];

% What kind of experiment did you do?
tags = {'put_one_tag_here','put_another_one_here'};
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
entry.datalocations = {'replace_this_with_the_entire_path_for_your_data'};
    % example: entry.datalocations =
    % {'/home/rebecca/matlabdata/experiemntaldata/2007_08_21/experiment1/'}
    %
    % If your data is in more than one place, you can include both of them:
    % entry.datalocations = {'/home/rebecca/foldername/date/exp1/',...
    % 'home/rebecca/foldername/date/exp2/'};
    
%% (blah blah - put options fields walk-through here)

%% Validate and save!

% put walk through of validation, and goign through the error_messages that
% you get, and how to fix common problems...

% save step: picking a name and saving

% mention where to go when you want to acutally search for something...