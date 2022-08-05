function [entry, xdb_id] = seea_xdb_gui(varargin)
% SEEA_XDB_GUI interfaces with 'seea' suite to assist making xdb file
%
% Syntax: [entry, xdb_id] = seea_xdb_gui(varargin)
%   By itself, seea_xdb_gui will return 'entry' using the current directory
%        (if possible)
%  'varargin' supplied as struct 'entry' followed by optional paired arguments (for now):
%    -> 'directory' - full path: points to a directory with single electrode ephys data
%    -> 'xdb_id' - string (user_number) xdb_id (either to an existing xdb_id or 'new'
%    -> 'analysis_file' - a string pointing to a valid '.analysis' file
% 
% This function will attempt to load any specific information about the
% file from the inputs given ('directory' or '.analysis' file for now) and
% will then load up a GUI interface with options for this experiment.
% Note: for identical experiments on a given day, multiple directories can
% be linked under a single 'xdb' entry.
%
% See also SINGLE_ELECTRODE_EPHYS_ANALYZE, VALIDATE_XDB_ENTRY_VAR

% Copyright 2008 Julian P. Meeks
%
% Version History:
% 2008_05_16--> in progress: (JPM)

%% Check varargin so we know what we're up against:
options = struct;    % structure to contain options set by varargin

if size(varargin,2) < 1
    fprintf('Will attempt to load current directory: %s\n',pwd);
    options.directory = pwd;
elseif size(varargin,2) == 1
    if isstruct(varargin{1})
        if validate_xdb_entry(entry)
            options.entry = varargin{1};
        else
            warning('Supplied entry does not pass validation!');
            options.entry = varargin{1};
        end
    else
        error('Single-argument calls must contain an ''entry'' structure variable.');
    end
else
    if isstruct(varargin{1})
        if validate_xdb_entry(entry)
            options.entry = varargin{1};
            start = 2;
        else
            warning('Supplied entry does not pass validation!');
            options.entry = varargin{1};
            start = 2;
        end
    else
        start = 1;
    end
    for i = start:2:size(varargin,2)
        if ~ischar(varargin{i})
            error('Paired arguments must start with a valid string.');
        end
        switch varargin{i}
            case 'directory'
                if ischar(varargin{i+1})
                    if exist(varargin{i+1},'dir')
                        options.directory = varargin{i+1};
                    else
                        error('%s is not a valid directory',varargin{i+1});
                    end
                end
            case 'xdb_id'
                if ~isstr(varargin{i+1})
                    error('Supplied .xdb entry must be a string.');
                elseif strmatch('new', varargin{i+1}, 'exact')
                    options.xdb_id = 'new';
                else
                    xdb_check = varargin{i+1};
                    uscore = strfind(xdb_check, '_');
                    username = xdb_check(1:uscore-1);
                    if xdb_check(end-3:end) == '.xdb'
                        xdb_check = xdb_check(1:end-4);
                    end
                    fullpath_xdb_check = ['/usr/lab/exp_db/' username filesep xdb_check '.xdb'];
                    if ~exist(fullpath_xdb_check,'file')
                        error('%s is not a valid ''.xdb'' file',fullpath_xdb_check);
                    else
                        options.xdb_id = xdb_check;
                    end
                end
            case 'analysis_file'
                if ~isstr(varargin{i+1})
                    if ~isstruct(varargin{i+1})
                        error('You must supply a string or analysis struct for ''analysis file''.');
                    else
                        options.analysis = varargin{i+1};
                    end
                else
                    analysis_check = varargin{i+1};
                    if ~exist(analysis_check,'file')
                        error('%s is not a valid filename',analysis_check);
                    else
                        load(analysis_check,'-mat');
                        if ~exist('analysis', 'var')
                            error('Invalid ''analysis'' file: must contain ''analysis'' structure.');
                        else
                            options.analysis = analysis;
                            clear analysis;
                        end
                    end
                end
        end
    end
end
                
options = default(options, 'directory', pwd);
options = default(options, 'entry', struct);
options = default(options, 'xdb_id', 'new');
options = default(options, 'analysis', struct);

%% Okay, no more of that 'boring' stuff, now onto setting up GUI fields!
valid_xdb = validate_xdb_entry_var;






end