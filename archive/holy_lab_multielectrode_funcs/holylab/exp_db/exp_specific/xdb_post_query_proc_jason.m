function entry=xdb_post_query_proc_jason(entry) 
   if(isequal(entry.tags, {'singing'}))
      oldStr='$EXP_DATA_ROOT$';
      newStr='@label=jason_behavior_e;path=jason/exp_data_4219to4624.full@';
      % or: newStr='@host=envy;path=/mnt/jason_behavior_e/jason/exp_data_4219to4624.full@';
      
      entry.data_locations=strrep(entry.data_locations, oldStr, newStr);
   elseif(isequal(entry.investigator, {'francesco'}))
      return
   else
      error('coding error?'); 
   end
 
   