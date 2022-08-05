% Script to work out converting Rebecca's dcell-style information into the
% general database format...

% Plan: 
%   - one entry per sorted cell, since that seems to be the fundamental
%   unit I'm working in these days...
%   - 


%data = load('/rocky/rebecca/data/pooledData/aob_small_2007_02_12.cells','-mat');

db_entries = struct;
for c = 1:length(data.aob)
    tc = data.aob(c);
%% add all the fields!
    % version
    db_entries(c).version = 1;
    % investigator
    db_entries(c).investigator = {'rebecca'};
    % start_date
    ds = tc.recordkeeping.date;
    db_entries(c).start_date = [str2num(ds(1:4)) str2num(ds(6:7)) str2num(ds(9:10))];
    % tags (try to put as many as are relevant...) --> still need to screen
    % for lfp, heartrate!
    clear tags
    tags = {'aob','single_electrode','sorted','mitral_cell'};
    drugs_map = {'bic','bicuculline','bic?','bicuculline',...
        'dcg-iv','dcg-iv','DCG','dcg-iv',...
        'ttx','ttx','TTX','ttx','ttx?','ttx','TXX','ttx',...
        'RC','negative_control_drug'};
    if isfield(tc,'drug')
        if isfield(tc.drug,'used')
            if tc.drug.used
                if isfield(tc.drug,'ID')
                    ID = {tc.drug.ID};
                else % stupid specialized case
                    ID = {tc.drug.drug1.ID,tc.drug.drug2.ID,tc.drug.drug3.ID};
                end
                for iid = 1:length(ID)
                tags{end+1} = apply_string_map(ID{iid},drugs_map);
                if strcmp(tags{end},'no match found')
                    error('unrecognized drug type?')
                end
                end
            end
        end
    end
    db_entries(c).tags = tags;
    % stimuli
    utags = tc.uniquetags;
    nStim = length(utags);
    stim = cell(1,nStim);
    substance_map = {'R','negative_control','U','urine_whole','K','K'};
    conc_unit_map = {'R','none','K','molar','U','fold'};
    for s = 1:nStim
        clear tstim
        tstim = struct;
        tstim.user_tag = utags{s};
        td = tc.tagdefs{s}.def;
        tstim.duration = td.time;
        tstim.category = apply_string_map(td.substance,substance_map);
        tstim.identity = '';
        tstim.concentration = td.concentration;
        tstim.strain = td.strain;
        tstim.conc_unit = apply_string_map(td.substance,conc_unit_map);
        if strcmp(td.sex,'M+F')
            tstim(2) = tstim;
            tstim(1).sex = 'M';
            tstim(2).sex = 'F';
        else
            tstim.sex = td.sex;
        end
        stim{s} = tstim;
    end
    db_entries(c).stimulus = stim;
    % data_locations
    db_entries(c).data_locations = tc.path;
    % comment
    db_entries(c).comment = tc.analysisnotes.general;
%% validate and save...
    [valid error_messages] = validate_xdb_entry_var(db_entries(c));
    if ~valid
        stoppoint = 1;
    end
end

