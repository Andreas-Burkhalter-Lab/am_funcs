entries=search_xdb(@(entry) isfield(entry, 'investigator') && isequal(entry.investigator, {'francesco'}));

fns=cell(size(entries));
[fns{:}]=deal('id');

ids=cellfun(@getfield, entries, fns, 'uniformoutput', false);

cellfun(@delete_xdb_entry, ids);
