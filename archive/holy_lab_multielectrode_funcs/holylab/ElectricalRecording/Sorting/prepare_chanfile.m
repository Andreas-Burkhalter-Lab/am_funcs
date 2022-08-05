function chanfile = prepare_chanfile(dirname)
% PREPARE_CHANFILE: generate list of sorting decisions automatically
% Syntax:
%   chanfile = prepare_chanfile(sorting_dirname)
% where
%   sorting_dirname is the name of the CASS-sorted directory
% and
%   chanfile is a structure array of channel#/sortingfile pairs, i.e., an
%     array of structures with fields "channel" and "file".  This is the
%     format expected by EPHYS/EPHYSFETCH. By default, this looks for
%     files of the name "1.mat" or "3.mat", and chooses the one with the
%     highest number.  It also warns you when there is any
%     ambiguity. Directories that do not contain sorting files with these
%     names are omitted from the chanfile list.
%
% See also: SORTED_GET_CHANNELS.
  
% Copyright 2007 by Timothy E. Holy
  
  dcur = pwd;
  cd(dirname)
  % Parse channel #s
  chandirs = dirbyname('chan*');
  for i = 1:length(chandirs)
    channum(i) = str2num(chandirs{i}(5:end));
  end
  [channum,sortIndex] = sort(channum);
  chandirs = chandirs(sortIndex);
  % Successively traverse directories
  dtop = pwd;
  chanfile = [];
  for i = 1:length(chandirs)
    cd(chandirs{i})
    filenames = dirbyname('*.mat');
    validnames = {};
    validnums = [];
    % Match the names that begin with a numeric tag
    for fileIndex = 1:length(filenames)
      tmp = sscanf(filenames{fileIndex},'%d',1); % does it start with a #?
      if ~ischar(tmp)
        validnames{end+1} = filenames{fileIndex}(1:end-4); % clip extension
        validnums(end+1) = tmp;
      end
    end
    if ~isempty(validnums)
      [mxnum,mxIndex] = max(validnums);
      if (length(validnums) > 1)
        warning('chanfile:ambiguous',['On ' chandirs{i} ...
          ', multiple sorting files were found. Choosing ' ...
          validnames{mxIndex} '.mat.']);
      end
      chanfile(end+1).file = validnames{mxIndex};
      chanfile(end).channel = channum(i);
    else
      warning('chanfile:none',['No sorting results found on ' chandirs{i}]);
    end
    cd(dtop)
  end
  cd(dcur)
  
