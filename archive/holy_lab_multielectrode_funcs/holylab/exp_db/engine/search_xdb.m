function result=search_xdb(callback, options)
% SEARCH_XDB: search the experiment database
% Syntax:
%   result=search_xdb(callback)
%   result=search_xdb(callback, options)
%   
% PRE:
%    callback: a function handle whose input arg is an entry and whose
%              output is true when the entry is desired.
%              SEE: query_callback_behavior.m
%    options: current only one field is supported:
%       .dbroot='/usr/lab/exp_db': the root dir of the database
% POST:
%    result: a cell array of database entries
%
% See also: QUERY_CALLBACK_BEHAVIOR.

   default_options('dbroot', '/usr/lab/exp_db');
   
   ownerDirs=dir(options.dbroot);
   
   result={};
   for idxOwner=1:length(ownerDirs)
      owner=ownerDirs(idxOwner).name;
      if(isequal(owner, '.') || isequal(owner, '..'))
         continue; 
      end
      
      entryFiles=find_first_pattern({'*.xdb'}, fullfile(options.dbroot, owner));
      if(isempty(entryFiles))
         continue; 
      end
      
      for idxEntry=1:length(entryFiles)
         entry=load(entryFiles{idxEntry}, '-mat');
         
         % matched=false;
         try
            if(iscell(callback))
               matched=callback{1}(entry, callback{2:end});
            else
               matched=callback(entry);
            end
         catch
            matched=false;
         end
         
         if(matched)
           try
            ownerPQP=str2func(['xdb_post_query_proc_' owner]);
            entry=ownerPQP(entry);
           end
           result{end+1}=entry;
         end % if, qualified
      end
   end % for, each owner
   
   
