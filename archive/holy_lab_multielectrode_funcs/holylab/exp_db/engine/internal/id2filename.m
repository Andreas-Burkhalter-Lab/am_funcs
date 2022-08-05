function filename=id2filename(id, options)
   default_options('dbroot', '/usr/lab/exp_db');
   
   owner=id2owner(id);
   
   filename=fullfile(options.dbroot, owner, [id '.xdb']);
   
