function xdb_ids = get_xdb_ids_from_entry(entry)
% GET_XDB_IDS_FROM_ENTRY
% 
% function to parse through a linear cell array of entries (structs with
% potentially non-overlapping field names), returning only the unique ids
% from each entry
%
% SEE ALSO: GETFIELD
%
if isfield(entry, 'id')
    xdb_ids = entry.id;
else
    xdb_ids = NaN
end