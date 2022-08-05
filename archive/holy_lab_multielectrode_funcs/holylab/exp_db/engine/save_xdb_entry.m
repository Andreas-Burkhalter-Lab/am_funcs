function id=save_xdb_entry(entry, options)
% save experiment into the db.
% SYNTAX:
%    save_xdb_entry(entry, options)
% PRE:
%    entry: the data to save
%    options: an optional struct with fields: dbroot, id, owner, retries
% NOTE: 
%    if you pass both id and owner, owner is discarded.
   
   default_options('dbroot', '/usr/lab/exp_db');
   default_options('id', []);
   default_options('owner', '');
   default_options('retries', 2); % # of retries when can't lock the entry
   
   if(isempty(options.id))
      if(isempty(options.owner))
         owner=getenv('USER');
      else
         owner=options.owner;
      end
   else
      owner=id2owner(options.id);
   end
   
   if(isempty(options.id))
      options.id=gen_xdb_entry_id(owner, struct('dbroot', options.dbroot));
   end
   
   entry.id=options.id; % NOTE: the old id, if exists, will be replaced
   
   filename=id2filename(entry.id);
   
   if(~trylock_file(filename, struct('retries', options.retries)))
      error(['can''t acquire the lock to the database entry: ' entry.id]);
   end
   
   save(filename, '-struct', 'entry', '-mat');
   
   unlock_file(filename);
   
   id=entry.id;
   
   
      