function analysis = single_electrode_ephys_analyze(varargin)
% SINGLE_ELECTRODE_EPHYS_ANALYZE analyzes single electrode data from directory
% 
% Note: 'seea( )' is a function shorthand for single_electrode_ephys_analyze
% Usage: analysis = seea( ) -- analyzes current directory only
%        analysis = seea(...,'directory', directory) -- analyzes directory indicated
%        analysis = seea(...,'recursive', true/false) -- recurses through directories
%        analysis = seea(..., 'manual_valvelabel', _filename_)
%                   -- allows use of input filename (function) to assign valves.  This
%                      is important on recordings where errors in automatic 
%                      valve assignment are evident "filename" must refer
%                      to a function which returns a cell array of strings
%                      including the desired valve labels.
%        analysis = seea(...,'call_func', @function)
%                   -- allows one to designate a single function or "function of functions"
%                      that take as inputs an ephys structure and returns
%                      an updated "analysis" structure with fields related
%                      to the analysis set
%                      By default, seaa calls '@seea_deltar' which calculates
%                      some basic properties of the data set, namely firing
%                      rates, changes in firing rate, peak firing rates,
%                      etc.
%                      example seaa(...,'call_func',
%                      @seea_lifetime_sparseness) must return 
%                      {struct}{struct}...{struct} up to size(ephys,2)
%        analysis = seea(...,'defined_subintervals', {'s%de%d'})
%                      allows user to flag analysis structure such that 
%                      subsequent calls by analysis functions will perform
%                      their analyses on specific subintervals within the
%                      analysis structure
%        analysis = seea(...,'exp_db', true/false)
%                   -- (in progress) allows 
%                   
% 
% Precursor to a more robust analysis package for ephys data for the
% sulfated steroid tuning project.  This script/function should eventually:
%   <A> load data from the current directory (or supplied directory with
%       potential for recursive seeking)
%   <B> identify all unique compounds and concentrations for stimuli
%   <C> Do as extensive an analysis as we can imagine, saving data in a
%       "analysis.mat" file in the raw data directory.  Include:
%          <a> delta r values using sliding window and various integration times
%          <b> selectivity values using sliding window/various integration
%          <c> discriminability values using sliding window/ etc.
%          <d> lifetime sparseness values ...
%          <e?> direct analysis of drug intervals
%          <f> t-test comparisons for each interval compared with controls
%          <g> ranksum comparisons for each interval...
%          <h> a per-cell "tag" to identify as putative/confirmed cell type
%              (mitral, periglomerular, granule, etc.)
%   <D?> Give user the potential to add items to the data for saving as an
%        .xdb file (strain, age, tags, stimuli, tags, etc.)
%
%   Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)

% Revision history:
% 2008_04_28-2008_05_09: Wrote it (JPM)

%% Parse inputs (varargin)

options = struct;  % structure variable containing the function options

if size(varargin{1},2) < 2
    if size(varargin{1},2) == 1
        errormsg = 'Please use paired input arguments.  See single_electrode_ephys_analyze for usage.';
        error(errormsg);
    end
else
    for arg_idx = 1:2:size(varargin{1},2)
        if ~ischar(varargin{1}{arg_idx})
            errormsg = 'Odd inputs must be strings.  See single_electrode_ephys_analyze for options.';
            error(errormsg);
        else switch varargin{1}{arg_idx} 
            case 'directory'
                if ischar(varargin{1}{arg_idx+1})
                    if exist(varargin{1}{arg_idx+1},'dir')
                        options = default(options,'directory', varargin{1}{arg_idx+1});
                    else
                        errormsg = 'The input directory does not exist.  Please check typing';
                        error(errormsg);
                    end
                else
                   errormsg = 'Input directory must be a string.';
                   error(errormsg);
                end
            case 'recursive'
                if islogical(varargin{1}{arg_idx+1})
                    options = default(options,'recursive', varargin{1}{arg_idx+1})
                else
                    errormsg = '''recursive'' variable must be a logical';
                    error(errormsg);
                end
            case 'manual_valvelabel'
                if ~iscell(varargin{1}{arg_idx+1})
                    errormsg = 'The input ''manual_valvelabel'' must return a cell array of strings.'
                    error(errormsg);
                else
                    options = default(options, 'manual_valvelabel', varargin{1}{arg_idx+1});
                end
            case 'call_func'
                options.call_func = varargin{1}{arg_idx+1};
            case 'defined_subintervals'
                if iscell(varargin{1}{arg_idx+1})
                    for idx = 1:size(varargin{1}{arg_idx+1},2)
                        if ~ischar(varargin{1}{arg_idx+1}{idx})
                            errormsg = 'You must supply ''defined_subintervals'' as a cell array of strings or char array.'
                            error(erromsg);
                        elseif isempty(strfind(varargin{1}{arg_idx+1}{idx},'s'))||...
                               isempty(strfind(varargin{1}{arg_idx+1}{idx},'e'))
                           errormsg = 'You must supply ''defined_subintervals'' using s_e_ notation (start/end)';
                           error(errormsg);
                        end
                    end
                elseif ~ischar(varargin{1}{arg_idx+1})
                    errormsg = 'You must supply ''defined_subintervals'' as a cell array of strings or char array.'
                    error(erromsg);
                end
                options.defined_subintervals = varargin{1}{arg_idx+1};
            end
        end
    end
    options = default(options, 'directory', pwd);        % default is current directory (PWD)
    options = default(options, 'recursive', false);      % default recursive is FALSE
    options = default(options, 'manual_valvelabel', []); % default valvelabel set to EMPTY
    options = default(options, 'call_func', @seea_deltar);         % default call_func set to EMPTY
    options = default(options, 'defined_subintervals', 'all');   % default defined_subintervals chooses all
end

% Assert: options.directory
%                .recursive
%                .manual_valvelabel
% are set!

%% Load data from directory into ephys (from generic epanalyze1.m)
% gather names of the utilized snippet files in the current directory
% (*clean.ssnp or *.ssnp)
if ~isempty(dirbyname([options.directory filesep '*clean.ssnp']))
    files = dirbyname([options.directory filesep '*clean.ssnp']);
elseif ~isempty(dirbyname([options.directory filesep '*.ssnp']))
    files = dirbyname([options.directory filesep '*.ssnp']);
else
    errormsg = fprintf('There are no *.ssnp files in %s directory. Please snippet first.', options.directory);
    error(errormsg);
end
% for f_idx = 1:size(files,2)
%     files{f_idx} = [options.directory filesep files{f_idx}];
% end


% gather names of the .vlv files from the directory
if ~isempty(dirbyname([options.directory filesep '*.vlv']))
    stimfilenames = dirbyname([options.directory filesep '*.vlv']);
else
    errormsg = fprintf('There are no *.vlv files in %s directory.', options.directory);
    error(errormsg);
end
    
% in case naming convention does not load files in time order, force the issue:
%    files = sort_ai_by_time(files);

% initialize the ephys structure for the files in this directory 
    data = ephysfromai(files, struct('usefullpath',true));

% update data with the stimulus information (timing, valve numbers)
    data = ephysfetch(data,'stimulus');
    data = ephysfetch(data,'sniptimes');
    
% create a cell array of the available valve labels from the ephys
% structure:
%    data.valvelabels(1) = {'50 mM KCl'}; data.valvelabels(16) = []; % this expt had an error in the user comment section
%    [data.valvelabels] = deal(jpm_simple_mixed_urine_valvelabels_2007_04_01);  % UNCOMMENT FOR MANUAL VALVE LABELING
%    data.valvelabels(1) = {'Ringer''s'};  % mistake on 2008_03_06

% if exist([options.directory filesep 'analysis.mat'], 'file')
%     load('analysis.mat', '-mat');
%     if ~exist('analysis', 'var')
%         error('Existing ''analysis.mat'' file is not seea_compatible, consider changing name.');
%     else
%         [data.valvelabels] = deal(analysis.valvelabels);
%     end
% end

if ~isempty(options.manual_valvelabel)&&~exist('analysis', 'var')
    [data.valvelabels] = deal(options.manual_valvelabel);
end

if exist('manual_valvelabels.mat', 'file')
    load('manual_valvelabels.mat', '-mat');
    data.valvelabels = manual_valvelabels;
end

foundvalve = 0;
for valve_idx = size(data(1).valvelabels,2):-1:1
    if isempty(data(1).valvelabels{valve_idx})&&foundvalve == 0
        for data_idx = 1:size(data,2)
            data(data_idx).valvelabels(valve_idx) = [];
        end
    elseif isempty(data(1).valvelabels{valve_idx})&&foundvalve == 1
        for data_idx = 1:size(data,2)
            data(data_idx).valvelabels(valve_idx) = {'EMPTY'};
        end
    else
        foundvalve = 1;
    end
end
valvelabels = {data(1).valvelabels};
     
% set the names for the subdirectories containing sorted data information

[data.sortfile] = deal('sort1');
    chanfile.channel = 1; % the channel number to sort
    chanfile.file = '1'; % the basname of the .sorted file to load

% update data using 'ephys_cass_fetch'
    data = ephys_cass_fetch(data, chanfile);
    [data.sort_cass_chanfile] = deal(chanfile);
    data = ephysfetch(data, 'celltimes');
    
% 'trange' contains a 2-vector (in seconds) which delimits the time 
% before and after your stimulus onset to consider for plotting:
    trange = [-10 20];

% call 'intervalsfromstim' to assign arrays containing the intervals,
% valve identities, and valve names from this experiment
    [intervals,identities,vlvsout] = intervalsfromstim(data,trange);

% call ephyssubrange to update your data (ephys struct) to onlyfunction handle
% contain the sniptimes, etc from your designated intervals
    
    data = ephyssubrange(data,intervals);
    
% one-line 'if' statement to collapse your data across cells into a
% one-dimensional cell array:
    if iscell(data), data = cat(2,data{:}); end

% call 'deal' (MATLAB function) to add the offset time from trange
% to every cell in the data array (this allows the time of stimulus
% onset to be equal to zero)
    [data.toffset] = deal(trange(1));

% collapes the valve identities corresponding to each
% interval into a 1-dimensional cell array
    identities = cat(2,identities{:});

% uncomment the following lines to use stimulus timing in the tags
    %tagops.addtime = 1;         % flag to add the time
    %tagops.timedelim = ' $ ';   % will occur on plot before numeral
    %tagops.timeprec = 0.5;      % 0.5 sec precision (sloppy)
    
% set the "tag" of each interval to the name of the valve label
% 'utags' will become a cell array of strings from all used valves
    %[data,utags] = ephystag(data,valvelabels{1}(identities), tagops);
    [data,utags] = ephystag(data,valvelabels{1}(identities));
    
%% Prepare analysis structure for passing to functions
% Load in any existing 'analysis.mat' file in directory
if exist([options.directory filesep 'analysis.mat'],'file')
    load([options.directory filesep 'analysis.mat'],'-mat');
else
    analysis = struct;
end

% Add valvelabels to analysis structure:
if ~isfield(analysis,'valvelabels')
    [analysis.valvelabels, link, first] = unique(valvelabels{1}, 'first');
end

% Add directory as identifier
if ~isfield(analysis,'id')
    analysis.id = options.directory;
end

% Add/update 'defined_subintervals'
if ~iscell(options.defined_subintervals)
    if options.defined_subintervals == 'all'
        if ~isfield(analysis, 'defined_subintervals')
            analysis.defined_subintervals = options.defined_subintervals;
        end
    end
else
    analysis.defined_subintervals = options.defined_subintervals;
end

%% Call functions using ephysin and analysis
if iscell(options.call_func)
    for call_idx = 1:size(options.call_func,2)
        analysis = options.callf_func{call_idx}(data,analysis);
    end
else
    analysis = options.call_func(data,analysis);
end

%% save data into directory
save([options.directory filesep 'analysis.mat'], 'analysis');

% to-do: set up a system whereby a copy is saved with a systematic file
% identifier to /usr/lab/exp_db/ if desired

end