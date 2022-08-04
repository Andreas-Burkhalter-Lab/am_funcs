function varargout = saveIfExists(filename, varsIn, varargin)
%SAVEIFEXISTS 'Save.m', ignoring nonexistent variables. 
%%% Inputs:
%%%%% 1. filename = save to save variables into
%%%%% 2. varsIn = cell array of strings naming variables to save
%%%%% 3. varargin = other arguments to pass to 'save.m', such as version and
%%%%%       '-append'
%%% Optional output argument will list the existent variables from varsIn
%%%     which were saved. 
%%% last edited 20/22/15

varsInputToSave = cell(size(varsIn)); % variables to input to save.m

for indVar = 1:numel(varsIn)
    existExp = sprintf('exist(''%s'',''var'')',varsIn{indVar}); % to check for existence in base workspace
    if evalin('caller',existExp) % if it exists in caller workspace
        varsInputToSave{indVar} = varsIn{indVar}; % if variable exists, input to save.m
    end
end

varsInputToSave = varsInputToSave(~cellfun(@isempty,varsInputToSave)); % remove empty cells

if nargout == 1
    varargout{1} = varsInputToSave;
elseif nargout > 1
    error('Only one output argument may be specified.')
end

varsString = [];
for indVar = 1:numel(varsInputToSave)
    varsString = [varsString '''' varsInputToSave{indVar} ''','];
end
varsString = varsString(1:end-1); % remove last comma

% Make and execute save expression in caller workspace. 
if isempty(varargin)
    saveExpression = sprintf('save(''%s'',%s)',filename,varsString); % no extra optional arguments to save.m
else
    if iscell(varargin{1})
        varargin = varargin{:}; % in case varargin was passed as a cell of strings, rather than list of strings
    end
    argString = [];
    for indVar = 1:numel(varsInputToSave)
        argString = [argString '''' varargin{indVar} ''','];
    end
    argString = argString(1:end-1); % remove last comma
    saveExpression = sprintf('save(''%s'',%s,%s)',filename,varsString,argString);
end

evalin('caller',saveExpression); % execute 'save.m' on existent variables with optional parameters

