function get_xdb_entry_template(op)

% Opens template to walk user through adding entry to xdb database
%
% SYNTAX:   get_xdb_entry_template
%           get_xdb_entry_template(op)
%
% Just opens the template found at "template_for_adding_xdb_entry_source.m"
% and saves it as "template_for_adding_xdb_entry.m" in the current
% directory, so that it can be played with and modified as needed.
%
% OPTIONS:
%   open_option = 1 (default) or 0; just means opens the script as well as
%                   saving it in the current directory
%
% HISTORY:
%   2007-08-21  (RCH)   wrote it; based loosely off getscripts.m
%
% See also: template_for_adding_xdb_entry_source,
%   display_xdb_entry_options, validate_xdb_entry_var

%% input and options

if nargin < 1
    op = struct;
end
op = default(op,'open_option',1);

%% do it...

% figure out what script were dealing with
scriptname = 'template_for_adding_xdb_entry';
sourcename = [scriptname '_source'];
targetname = [scriptname '.m'];

% read in the info
hs = which(sourcename);
fid = fopen(hs);
content = fread(fid,'*char');
fclose(fid);

% now, write it
[fid, message] = fopen(targetname,'w');
fwrite(fid, content, 'char');
fclose(fid);
if op.open_option
    evalstring = ['open ' targetname];
    eval(evalstring)
end
