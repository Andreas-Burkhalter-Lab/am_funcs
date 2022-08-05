function display_xdb_entry_options(fieldname)

% Displays options for fields in a valid xdb database entry
%
% SYNTAX:   display_xdb_entry_options
%           display_xdb_entry_options(fieldname)
%        
% This function just displays formatted information on what fields are
% currently allowed in a database entry, and what restrictions there are on
% what can be put in those fields.  It is primarily a wrapper for
% validate_xdb_entry_var.  If a particular fieldname is passed, it will
% give information on that field in particular; otherwise, it will just
% tell you everything it knows about how to make an entry valid...
%
% HISTORY:
%   2007-08-21  (RCH)   wrote it
%
% See also: validate_xdb_entry_var, get_xdb_entry_template

%% options

if nargin < 1
    do_all = 1;
    fieldname = [];
else
    do_all = 0;
end

%% get current info from validate_xdb_entry_var

valid = validate_xdb_entry_var;

%% format the results and print to command line for the user

if do_all

    fprintf(1,'Required fields:\n')
    
    fprintf(1,'oops - haven''t finished writing this, until format settles down a bit...\n')
    
else
    
    fprintf(1,['\nThe field "' fieldname '" can have\none or more entries from the following list:\n\n     '])
    entry_options = valid.limited_field_options.(fieldname);
    nEntries = length(entry_options);
    max_per_line = 4;
    for e = 1:nEntries
        fprintf(1,entry_options{e})
        if e<nEntries
            fprintf(1,', ')
        end
        if (e == nEntries)
            fprintf(1,'\n')
        elseif (mod(e,4) == 0)
            fprintf(1,'\n     ')
        end
    end
    fprintf(1,['\nIf something needs to be added to this list,\n'...
        'it can be done by opening validate_xdb_entry_var.m\n'...
        'and editing the top section of code\n\n']);
    
end
