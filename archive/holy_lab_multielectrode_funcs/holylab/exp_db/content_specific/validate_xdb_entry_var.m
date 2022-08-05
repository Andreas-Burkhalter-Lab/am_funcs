function [valid error_messages] = validate_xdb_entry_var(entry)

% this function checks if an experiment database entry variable is valid
%
% SYNTAX: [valid error_messages] = validate_xdb_entry_var(entry)
%         [requirements_for_validity] = validate_xdb_entry_var
%
% If no entry is passed to the function, it will return information on what
% requirements there are for validity
%
% HISTORY:
%   2007-08-??  (ZG)    wrote it
%   2007-08-20  (RCH)   added restrictions for certain fields with limited
%                       options, and someoutput to help ID errors
%   2007-08-21  (RCH)   added in structural changes suggested at lab
%                       meeting
%   2007-08-27  (JPM)   added additional tags for experimental preparation
%                       and histology/imaging for fluoro_jade
%   2007-10-5   (TEH&JPM) debugging & improvement in error messages

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SETTINGS: The basic structure of a database entry.  This is where you
%%% can both read about what fields an entry has, as well as change what
%%% entries are allowed...


%% What are the required and optional fields?

required_fields = {...
    'version','investigator','start_date','tags',...
    'data_locations'...
    };
optional_fields = {...
    'stimulus','comment','analysis',...
    'experimental_animal','related_entries'...
    'section', 'histology', 'vno_aob_hemihead',...
    'confocal_imaging', 'fluoro_jade',...
    'hist_analysis','id',...
    };

    % note that as much as possible, "category" divisions should be placed
    % in two existing tags fields, and new optional_fields should be reserved
    % for things that need more complex structure or named context
    % information....
    
    
%% "tags" options

% This field may contain as many entries as you'd like off the following
% list...
limited_field_options.tags = {...
     ... % I. basic experiment types - should probably make sure you have at least one of these!
        'single_electrode','mea','behavior','histology','imaging',...
     ...    % refinements on these - probably want both the general and specific terms, here
     ...     % 1) electrophys refinements
                    'lfp','sorted','unsorted','heart_rate',...
     ...     % 2) behavior refinements
                    'singing','bruce',...
     ...     % 3) imaging refinements
                    'confocal_imaging','fluoro_jade','section', 'hist_analysis'...
     ...     % 4) experimental preparation refinements
                    ,'vno_aob_hemihead','anesthetized',...
     ... % II. where/what for imaging or electrophys, are you recording from?
            'aob','vno','mitral_cell','granule_cell',...
     ... % III. other info that might help someone, eg drugs that might have been used
            'bicuculline','dcg-iv','ttx','negative_control_drug'};
% (if something new needs to be added, adding it here will carry out all
% necessary changes for it to be allowed by the validation function and
% listed by the display_xdb_entry_options function)
 

%% "investigator" options
        
% This field may contain as many entries as you'd like off the following
% list.  The first one is presumed to be the "primary" experimentor.
limited_field_options.investigator = {'jason','rebecca','francesco','diwakar',...
    'sandra','james','terry','xiaoyan','tim','julian', 'hannah'};


%% mouse_strain options (for both "stimulus" and "experimental_animal" fields)

mouse_strains = {'Balb/c','cba','B6D2F1/J'};


%% "experimental_animal" options

% I.  what fields are allowed at all?
limited_subfields.experimental_animal.fields = {...
    'age','sex','strain','tags'
    };

% II. what are the field-specific restrictions?
% 1) sex -
limited_subfields.experimental_animal.field_options.sex = {...
    'M','F'};
% 2) strain - (these were set above, under "mouse_strain options")
limited_subfields.experimental_animal.field_options.strain = ...
    mouse_strains;
% 3) tags -
limited_subfields.experimental_animal.field_options.tags = {...
    'TRP-/-','M','F',...
    'socially_experienced','socially_naieve',...
    'castrated'...
    };
% 4) age -
limited_subfields.experimental_animal.field_formats.age = ...
    'numeric';

%% "stimulus" options

% I. what fields are allowed?
limited_subfields.stimulus.fields = {...
    'category','identity','concentration',...
    'conc_unit','duration','user_tag',...
    'strain','sex', 'procedure'};
    % want to make sure have the right fields to begin with!
    
% II. what are the field-specific restrictions?
% 1) category -
limited_subfields.stimulus.field_options.category = {...
    'urine_whole','urine_fraction','K','negative_control','synthetic_compound'};
    % note: this field should probably stay fairly restrictive...
% 2) conc_unit -
limited_subfields.stimulus.field_options.conc_unit = {...
    'relative','mM', 'pM', 'nM', '\muM','none'}; % ('fold' means dilution factor)
% 3) strain - (these were set above, under "mouse_strain options")
limited_subfields.stimulus.field_options.strain = ...
    mouse_strains;
% 4) sex -
limited_subfields.stimulus.field_options.sex = ...
    {'M','F'};
% 5) concentration -
limited_subfields.stimulus.field_formats.concentration = ...
    'numeric';
% 6) identity -
limited_subfields.stimulus.field_formats.identity = ...
    'string';
% 7) duration -
limited_subfields.stimulus.field_formats.duration = ...
    'numeric'; % sould be in SECONDS
% 8) user_tag -
limited_subfields.stimulus.field_formats.user_tag = ...
    'string';


%% "analysis_record" options

% I. what fields are allowed?
limited_subfields.analysis.fields = {...
    'tags','location/filename','analysis_options',...
    'save_options','date'};

% II. what are the field-specific restrictions?
% (none yet)

%% "related_entries" options

% I. what fields are allowed?
limited_subfields.related_entries.fields = {...
    'UID',... % this refers to the uniqe identifier of the entry, ie, the name it's saved under
    'tags'};  % these tags are indexed to match the UID's, and can provide referential information - eg, things like 'parent'

% II. what are the field-specific restrictions?
% (note yet)

%%%%% END SETTINGS SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%% if there's no entry...
%           ... just return information on what it takes to be a valid entry

if nargin == 0
    error_messages = cell(1,0);
    valid.required_fields = required_fields;
    valid.optional_fields = optional_fields;
    valid.limited_field_options = limited_field_options;
    valid.limited_subfields = limited_subfields;
    return
end
    
%% do the checks...

    error_messages = cell(1,0);

    % 1) first, is the basic structure correct?
    
    fieldsused = fieldnames(entry);
    % are all the required fields there?
    
    if any(~check_strings_on_list(required_fields,fieldsused))
        i_missingfield = find(~check_strings_on_list(required_fields,fieldsused));
        error_messages{end+1} = ['You are missing one or more required fields: ' ...
            required_fields{i_missingfield}];
    end
    % are all the fields there at least options?
    if any(~check_strings_on_list(fieldsused,{required_fields{:},optional_fields{:}}))
        i_newfield = find(~check_strings_on_list(fieldsused,[required_fields,optional_fields]));
        errmsg_tmp = 'You have (an) unrecognized field(s): ';
        errmsg_tmp = [errmsg_tmp sprintf('%s ',fieldsused{i_newfield})];
        error_messages{end+1} = errmsg_tmp;
    end
    
    % 2) do base fields with restrictions on their formats meet the
    % requirements?
    
    limited_fields = fieldnames(limited_field_options);
    for f = 1:length(limited_fields)
        this_entrys_field = entry.(limited_fields{f});
        % make sure it's a cell, so can give informative error if it's not
        if ~iscell(this_entrys_field)
            error_messages{end+1} = ['field "' limited_fields{f} '" needs to be a cell array...'];
        end
        if any(~(check_strings_on_list(this_entrys_field,limited_field_options.(limited_fields{f}))));
            i_unrecognized_entry = find(~(check_strings_on_list(this_entrys_field,limited_field_options.(limited_fields{f}))));
            error_messages{end+1} = ['field "' limited_fields{f} '" has one or more unrecognized entries: ' ...
                this_entrys_field{i_unrecognized_entry} '; if a new entry type needs to be added, please' ...
                ' add it to the approrpiate list at the beginning of validate_xdb_entry_var'];
        end
    end
    
    % 3) do subfields with restrictions on their formats meet the
    % requirements?
    
    lim_subfield_bases = fieldnames(limited_subfields);
    for f = 1:length(lim_subfield_bases)
        if isfield(entry,lim_subfield_bases{f})
            base_req = limited_subfields.(lim_subfield_bases{f});
            as_used = entry.(lim_subfield_bases{f});
            if ~iscell(as_used)
                as_used = {as_used}; % since some are in cells, do all cell-by-cell
            end
            % see if all fields used match allowed fields
            fields_allowed = base_req.fields;
            for nthEntry = 1:length(as_used)
                fields_used = fieldnames(as_used{nthEntry});
                if any(~(check_strings_on_list(fields_used,fields_allowed)))
                    i_unrecognized_entry = find(~(check_strings_on_list(fields_used,fields_allowed)));
                    error_messages{end+1} = ['field "' lim_subfield_bases{f} '" has one or more unrecognized subfields: ' ...
                        fields_used{i_unrecognized_entry} '; if a new entry type needs to be added, please' ...
                        ' add it to the approrpiate list at the beginning of validate_xdb_entry_var'];
                end
            end
            % do the requried smaller checks one by one; do ones with
            % allowed string lists first
            if isfield(base_req,'field_options')
                subfields_with_limits = fieldnames(base_req.field_options);
                for sf = 1:length(subfields_with_limits)
                    strings_allowed = base_req.field_options.(subfields_with_limits{sf});
                    for nthEntry = 1:length(as_used)
                        if isfield(as_used{nthEntry},subfields_with_limits{sf})
                            strings_used = as_used{nthEntry}.(subfields_with_limits{sf});
                            if any(~(check_strings_on_list(strings_used,strings_allowed)))
                                i_unrecognized_entry = find(~(check_strings_on_list(strings_used,strings_allowed)));
                                error_messages{end+1} = ['subfield "' lim_subfield_bases{f} '.'...
                                    subfields_with_limits{f} '" has one or more unrecognized subfields: ' ...
                                    strings_used{i_unrecognized_entry} '; if a new entry type needs to be added, please' ...
                                    ' add it to the approrpiate list at the beginning of validate_xdb_entry_var'];
                            end
                        end
                    end
                end
            end
            % now do the ones with type restrictions
            if isfield(base_req,'field_formats')
                subfields_with_types = fieldnames(base_req.field_formats);
                for sf = 1:length(subfields_with_types)
                    type = base_req.field_formats.(subfields_with_types{sf});
                    if strcmp(type,'string')
                        type = 'char'; % because gosh darn it, I just can't stand telling when to use which term...
                    end
                    for nthEntry = 1:length(as_used)
                        if isfield(as_used{nthEntry},subfields_with_types{sf})
                            entry_value = as_used{nthEntry}.(subfields_with_types{sf});
                            if length(entry_value)>0
                                pass = NaN;
                                teststring = ['pass = is' type '(entry_value);'];
                                eval(teststring)
                                if isnan(pass)
                                    error('something went wrong in the code to check the type')
                                elseif ~pass
                                    error_messages{end+1} = [...
                                        'subfield "' lim_subfield_bases{f} '.' ...
                                        subfields_with_types{sf} '" has a type that does not match '...
                                        'the required type ( ' type ')'];
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    % check to make sure the date field is correctly formated 
    % (for now, enforcing vector format, so that can be fast for datenum
    % handling, but still human readable)
    datefield = entry.start_date;
    if (length(datefield)~=3)
         error_messages{end+1} = 'there are not three numbers in your start_date field - this field needs to be in date-vector format; see help in the datenum function';
    elseif ~(datefield(1)>1980)
         error_messages{end+1} = 'the first entry in your start_date needs to be a year, with all 4 numbers written out';
    elseif ~isbetween(datefield(2),[.9 12.1])
        error_messages{end+1} = 'the second entry in your start_date needs to be a month';
    elseif ~isbetween(datefield(3),[.9 31.1])
        error_messages{end+1} = 'the third entry in your start_date needs to be a day';
    end
    %%% This code is only necessary if we go back to using datenum form...
%     if ~isnumeric(datefield)
%         error_messages{end+1} = 'be sure that your start_date is in datenum format!  (see the help on the datenum function for assistance)';
%         valid = 0;
%     else % turn it into a date vector to make sure seems reasonable
%         v = datevec(datefield);
%         year_range = [1990 2050];
%         if ~isbetween(v(1),year_range)
%             error_messages{end+1} = ['unless this experiment is from before ' num2str(year_range(1)) ' or after ' ...
%                 num2str(year_range(2)) ', your start_date is in the wrong format (see datenum function for help)'];
%             valid=0;
%         end
%     end

    if length(error_messages) > 0
        valid = 0;
    else
        valid = 1;
    end
        
    
    %% subfunctions
    
    function [l_pass] = check_strings_on_list(strings,list)
        if ~iscell(strings)
            strings = {strings};
        end
        l_pass = ones(1,length(strings));
        for s = 1:length(strings)
            if ~any(strcmp(strings{s},list))
                if length(strings{s})>0
                    l_pass(s) = 0;
                end
            end
        end