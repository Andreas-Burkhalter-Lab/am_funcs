rootdir='/mnt/francesco_b01/francesco/francesco_physiology';

paths=genpath(rootdir);
paths=cellstr(split_str(paths, pathsep));

entries={};

for dirIdx=1:length(paths)
   filepath=fullfile(paths{dirIdx}, 'analyze1.m');
   if(~fileexist(filepath)) continue; end

   if(isequal(filepath,'/mnt/francesco_b01/francesco/francesco_physiology/spikeanalysis/analyze1.m'))
      continue;
   end
   
   % TODO: may need pass options.date for specific exp
   entry=gen_mea_entry_francesco(filepath);
   
   [ok, msg]=validate_xdb_entry_var(entry);
   if(ok)
      entries{end+1}=entry;
   else
     if iscell(msg)
       for msgIndex = 1:length(msg)
         disp(msg{msgIndex});
       end
     else
      disp(msg);
     end
      error('please fix gen_mea_entry_francesco()');
   end
end

ids=cell(size(entries));
for idxEntry=1:length(entries)
   ids{idxEntry}=save_xdb_entry(entries{idxEntry}, struct('owner', 'francesco'));   
end

idfile=[datestr(now, 30) '.xdb_ids'];
save(idfile, 'ids', '-mat');
disp('done');
