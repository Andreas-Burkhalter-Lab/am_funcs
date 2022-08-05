function [matchedFiles, matchedIndex]=find_first_pattern(filenames, parentDir)
% find_first_pattern: find the first file pattern that has existent file(s).
% USAGE:
%    [matchedFiles, matchedIndex]=find_first_patten(filenames, parentDir)
% PRE:
%    filenames: a cell array of file name patterns
%    parentDir (optional): the parent dir
% POST:
%    matchedIndex: the index of the first pattern that has existing file(s).
%                  -1 when no match.
%    matchedFiles: a cell array of the matching file names.
%                  {} when no match.

   if(nargin==2)
      for idx=1:length(filenames)
         filenames{idx}=fullfile(parentDir, filenames{idx});
      end
   end
   matchedFiles={}; matchedIndex=-1;
   for idx=1:length(filenames)
      d=dir(filenames{idx});
      nFiles=length(d);
      if(nFiles)
         matchedFiles=cell(1, nFiles);
         % [matchedFiles{:}]=deal(d.name);
         for tt=1:nFiles
            matchedFiles{tt}=replace_filename(filenames{idx}, d(tt).name);
         end
         matchedIndex=idx;
         return;
      end % if, found existent
   end % for, each file pattern
   
   
