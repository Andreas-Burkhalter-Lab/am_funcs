function id=gen_xdb_entry_id(username, options)
% generate entry's id

   default_options('dbroot', '/usr/lab/exp_db');

   % may use uuid, but it seems overkilled.
   
   max_id=10^6;
   
   while(true)
      id=sprintf('%s_%06d', username, floor(max_id*rand));
      
      if(~fileexist(id2filename(id, options)))
         break;
      end
   end
   
