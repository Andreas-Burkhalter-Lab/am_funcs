function cluster_events_raw_batch(merecFilenames, options)
% cluster_events_raw_batch:
% PRE: 
%    merecFilenames: a string or cell array of strings
%    options: the options passed to cluster_events_raw()
%       NOTE: options.fileToSave and options.channels
%             are filled/overwritten by this function.

   if(nargin==1)
      options=struct;
   end
   
   if(~iscell(merecFilenames))
      merecFilenames={merecFilenames};
   end
   
   for fileIdx=1:length(merecFilenames)
      curMerec=merecFilenames{fileIdx};
      for part=1:2
         options.fileToSave=replace_extension(curMerec, ...
            sprintf('_part%d.raw_cluster', part));
         options.channels.what='half';
         options.channels.value=part;
         cluster_events_raw(curMerec, options);
      end
   end % for, each input merec file
   