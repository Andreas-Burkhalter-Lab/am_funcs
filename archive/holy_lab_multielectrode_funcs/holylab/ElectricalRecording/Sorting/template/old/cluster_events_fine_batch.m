function cluster_events_fine_batch(merecFilenames, options)

% cluster_events_fine_batch:
%
% SYNTAX: cluster_events_fine_batch(merecFilenames, options)
%         cluster_events_fine_batch(merecFilenames)
%         cluster_events_fine_batch(options)
%         cluster_events_fine_batch
%
% PRE: 
%    merecFilenames: a string or cell array of strings
%    options: the options passed to cluster_events_fine()
%
% See also: cluster_events_fine

   getMerec = 0;
   if(nargin==1)
       if isstruct(merecFilenames)
           options = merecFilenames;
           getMerec = 1;
       else
           options=struct;
       end
   end
   if nargin == 0 
       options = struct;
       getMerec = 1;
   end
   if getMerec
       potentialFiles = dirbytime('*_part1.raw_cluster');
       if length(potentialFiles)>1
           warning('You seem to have multiple part1.raw_cluster files... I''m picking one to use essentially at random!');
       end
       filenametemp = potentialFiles{end};
       merecFilenames = filenametemp((1:end-18));
   end
   
   if(~iscell(merecFilenames))
      merecFilenames={merecFilenames};
   end
   
   for fileIdx=1:length(merecFilenames)
      curMerec=merecFilenames{fileIdx};
      for part=1:2
         rawClusterFile=replace_extension(curMerec, ...
            sprintf('_part%d.raw_cluster', part));
        if options.resume_on
            fineClusterFile=replace_extension(rawClusterFile, '.fine_cluster');
            if exist(fineClusterFile,'file')
                fileToUse = fineClusterFile;
            else
                fileToUse = rawClusterFile;
            end
        else
            fileToUse = rawClusterFile;
        end
         cluster_events_fine(fileToUse,options);
      end
   end
