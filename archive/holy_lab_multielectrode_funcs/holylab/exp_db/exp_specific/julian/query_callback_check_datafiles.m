function matched=query_callback_check_datafiles(entry)
% this query callback function takes a current entry and examines all current .xdb files
% in your exp_db subdirectory for entries with overlapping data file locations to the current
% entry.
% The function returns true if an experiment contains overlapping data
% files, and false if it does not
%
% Copyright Julian P. Meeks 
% Created 08/27/07
%    version_history
%
% SEE ALSO: SEARCH_XDB

   matched=0;
   
   data_locations = load_data_locations_from_cookie(pwd);
   
   if(isempty(data_locations)); return; end
   
 %  temp_entry_locations = cell2mat(entry.data_locations(:)); % buggy if all data_locations are not the same length!!!
   
   for idx_in = 1:size(data_locations,2)
       for idx_entry = 1:size(entry.data_locations,2)
           if ~isempty(strmatch(data_locations{idx_in}, entry.data_locations{idx_entry}))
                 matched = 1;
                 return
           end
       end
   end
   matched = 0;
end