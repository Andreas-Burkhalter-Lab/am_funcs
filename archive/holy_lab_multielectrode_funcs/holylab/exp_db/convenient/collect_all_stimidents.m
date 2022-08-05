function stimident = collect_all_stimidents(entry)
% COLLECT_ALL_STIMIDENTS: return all unique stimulus identities
% Syntax:
%   stimident = collect_all_stimidents(entry)
% where
%   entry is a cell array of database entries (or a single entry)
% and
%   stimident is a cell array of strings, containing the list of unique
%     stimulus identities used across all experiments.
  
% Copyright 2007 by Timothy E. Holy
  
  if ~iscell(entry)
    entry = {entry};
  end
  n_entries = length(entry);
  stimident = {};
  for entryIndex = 1:n_entries
    stim = entry{entryIndex}.stimulus;
    n_stimuli = length(stim);
    for stimIndex = 1:n_stimuli
      stimident{end+1} = stim{stimIndex}.identity;
    end
  end
  stimident = unique(stimident);
  