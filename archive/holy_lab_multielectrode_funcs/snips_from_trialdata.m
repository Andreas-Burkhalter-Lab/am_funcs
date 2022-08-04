function [ outsnips ] = snips_from_trialdata(trialdata, snipdata, trialnumbers, unit_names, merec_for_trialwave)
%SNIPS_FROM_TRIALDATA Get SNIPS from specified units on specified trials
%described in the trialdata dataset from the snipdata struct.
%%% [ outsnips ] = waves_from_trialdata(trialdata, snipdata, trialnumbers, unit_names)      
%   trialdata = same struct from ss_getspikes/rf_getspikes
%   snipdata = struct output from sorthead_from_raw
%   trial_numbers = indices of trials of interest from trialdata (chronological) 
%   unit_names = unit names (not indices); vector or cell of numbers
%      and/or strings; numbers will be turned to strings as 'ch[double';
%      must match unit names in snipdata.unit_names

%%% last edited 10/22/15 on vivid

%%% not yet implemented:
%   merec_for_trialwave: if supplied, also get the waveform for the entire trial for
%       each channel (matching one in the .merec file) specified by a unit
%       name that is numeric or starts with 'ch'

%% Reformat unit names and get unit indices
unit_inds = inf(size(unit_names)); % indices from snipdata.unitnames
if isnumeric(unit_names) % make unit names into strings if they are numeric
    unit_names = cellfun(@(x)['ch',num2str(x)],num2cell(unit_names),'UniformOutput',false);
end

% If snipdata doesn't contain 'unit names' field, use 'channels' field
% instead and turn to strings starting with 'ch'. 
if ~isfield(snipdata,'unitnames')
    snipdata.unitnames = cellfun(@(x)['ch',num2str(x)],num2cell(snipdata.channels),'UniformOutput',false);
end

for i = 1:numel(unit_names)
    % Check that this unit name is contained in snipdata.unit_names
    if ~any(strcmp(unit_names{i},snipdata.unitnames))
        error(['Unit name ''' unit_names{i} ''' not found in snipdata.unitnames.'])
    end
    % Get unit index from snipdata.unitnames
    unit_inds(i) = find(strcmp(unit_names{i},snipdata.unitnames));
end

%% Get snips from units of interest on trials of interest.
% The observation name will preserve the original trial number.
outsnips = trialdata(trialnumbers,:); % trial data from specified trials

for uindind = 1:length(unit_inds) % uindind = index within unit_inds, not values of unit_inds
    unitind = unit_inds(uindind); % index within snipdata.unitnames
    outsnips.tempvar = cell(size(outsnips,1),1); % make new variable in outwaves, then rename
    outsnips.Properties.VarNames{strcmp(outsnips.Properties.VarNames,'tempvar')} = [unit_names{uindind} '_snips'];
    outsnips.tempvar = cell(size(outsnips,1),1); % make new variable in outwaves, then rename
    outsnips.Properties.VarNames{strcmp(outsnips.Properties.VarNames,'tempvar')} = [unit_names{uindind} '_sniptimes'];
    for trialind = 1:length(trialnumbers) % ind in trialnumbers, not in the original trialdata dataset
%         thistrial = trialnumbers(trialind); % trial number from original trialdata dataset
        current_snip_inds = outsnips.snip_ind{trialind,unitind}; % snip indices (col) in snipdata.snips{unitind}
        outsnips{trialind,[unit_names{uindind} '_snips']} = snipdata.snips{unitind}(:,current_snip_inds);
        outsnips{trialind,[unit_names{uindind} '_sniptimes']} = snipdata.sniptimes{unitind}(current_snip_inds); % scans
    end
end

outsnips.spikes = []; % remove unneeded variables
outsnips.snip_ind = []; % remove unneeded variables
