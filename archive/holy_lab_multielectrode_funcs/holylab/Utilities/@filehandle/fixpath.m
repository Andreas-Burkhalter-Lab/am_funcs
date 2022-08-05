function [fh,requiredUserInteraction] = fixpath(fh,pathsToTry)
% FILEHANDLE/FIXPATH: check for existence & correct the path if necessary
% Syntax:
%   fh = fixpath(fh)
%   [fh,requiredUserInteraction] = fixpath(fh,pathsToTry)
% where
%   fh is a filehandle
%   pathsToTry is a cell array of strings, each containing a possible
%     path to the file. The abspathstr of the filehandle will be tried;
%     any other places you want checked (e.g., the current directory)
%     need to be listed explicitly (e.g., pathsToTry = {'.'} or {pwd}).
% and on output, the abspathstr of the filehandle will be updated. If
% this had to ask the user for help in locating the file, the
% requiredUserInteraction flag will be true. In cases where neither
% automatic nor user means for identifying the file succeeded, on output
% fh will be -1.
  
% Copyright 2007 by Timothy E. Holy
  
  requiredUserInteraction = false;
  fullName = fullfile(fh.abspathstr,fh.filename);
  if exist(fullName,'file')
    return;
  end
  if (nargin < 2)
    pathsToTry = {};
  end
  for i = 1:length(pathsToTry)
    fullName = fullfile(pathsToTry{i},fh.filename);
    if exist(fullName,'file')
      fh.abspathstr = pathsToTry{i};
      return;
    end
  end
  % OK, failed to find it automatically, resort to user query
  requiredUserInteraction = true;
  [filename,pathstr] = uigetfile('*',...
				 ['Please help me find ' fh.abspathstr],...
				 fh.filename);
  if ~filename
    fh = -1;
    return
  end
  fh.filename = filename;
  fh.abspathstr = pathstr;
