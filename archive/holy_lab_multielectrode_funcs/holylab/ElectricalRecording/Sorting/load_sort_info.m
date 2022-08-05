function sort_info = load_sort_info(filename)
% LOAD_SORT_INFO: load spike sorting information for one channel
% Syntax:
%   sort_info = load_sort_info(filename)
% where
%   filename is the name of the file;
% and
%   sort_info is a structure with fields documented in SAVE_SORT_INFO.
%
% See also: SAVE_SORT_INFO, AUTOSORT, CASS.
  
% Copyright 2005 by Timothy E. Holy
  
  load(filename)
  % Originally, these were binary files with their own format...
