function delete_xdb_entry(id, options)
% delete an xdb record
% SYNTAX:
%   delete_xdb_entry(id, options)
% PRE:
%   id: the id of the entry
%   options: an optional struct with optional fields:
%     .dbroot=/usr/lab/exp_db: the root of the database
%     .retries=2: # of retries

   default_options('dbroot', '/usr/lab/exp_db');
   default_options('retries', 2);
   
   filename=id2filename(id, options);
   
   if(~trylock_file(filename, struct('retries', options.retries)))
      error(['can''t acquire the lock to the database entry: ' id]);
   end
   
   delete(filename);
   
   unlock_file(filename);
   
