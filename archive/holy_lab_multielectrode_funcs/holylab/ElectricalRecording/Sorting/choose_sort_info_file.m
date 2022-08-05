function filename = choose_sort_info_file(dirchan,force_autosort_info)
% CHOOSE_SORT_INFO_FILE: automatic & manual selection of sort_info file
%
% Syntax:
%   filename = choose_sort_info_file(dirchan)
%   filename = choose_sort_info_file(dirchan,force_autosort_info)
% where
%   dirchan is the name of the directory containing the sort_info files
%   force_autosort_info is a flag which, if true, causes autosort_info.mat
%     to be chosen (if present)
% and
%   filename is the name of the chosen file.
%
% HISTORY:
%
%   2005-08-09 - bug in user-selected filename assignment fixed (RCH)
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 2)
    force_autosort_info = 0;
  end
  sort_info_files = dirbyname([dirchan '*.mat']);
  % Remove the proj.mat file
  projIndex = strmatch('proj.mat',sort_info_files,'exact');
  sort_info_files(projIndex) = [];
  % Remove any time_amp*.mat files
  taIndex = strmatch('time_amp',sort_info_files);
  sort_info_files(taIndex) = [];
  if isempty(sort_info_files)
    error(['There are no sorting files in directory ' dirchan]);
  end
  if (length(sort_info_files) == 1)
    % If there's only one file, return that
    filename = sort_info_files{1};
  else
    autosortIndex = strmatch('autosort_info.mat',sort_info_files, ...
      'exact');
    if (force_autosort_info && ~isempty(autosortIndex))
      filename = 'autosort_info.mat';
      return
    end
    if (length(sort_info_files) == 2 && ~isempty(autosortIndex))
      % If there are two files, and one is autosort_info.mat, return the
      % other
      filename = sort_info_files{setdiff([1 2],autosortIndex)};
    else
      % Must ask the user to choose
      filenumber = listdlg('ListString',sort_info_files,...
        'SelectionMode','single',...
        'PromptString','Choose sorting file');
      filename = sort_info_files{filenumber};
    end
  end