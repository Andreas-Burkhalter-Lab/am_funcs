function clean_path()
% clean_path: exclude inappropriate directories from the path
% Syntax:
%   clean_path
%
% This excludes directories named things like ".svn" from the matlab
% path. It also will exclude any directory named "Old".
% On unix systems, the "getmatlabpath" function (called by aliases like
% mymatpath) takes care of this for you.
  
   excludeList={'CVS', '\.svn', '@', 'Old'}; 

   paths=path;
   paths=cellstr(split_str(paths, pathsep));
   
   for idxPath=1:length(paths)
      curPath=paths{idxPath};
      indices=regexp(curPath, excludeList);
      if(~isempty(cat(2,indices{:})))
         rmpath(curPath);
      end
   end
   