function cass_sort_apply_AM(dirname,options,snipfile_override)
% CASS_SORT_APPLY: execute sorting decisions for a directory tree
%%% AM edited 8/20/15 on vivid
%%% -use sortheader_importchan_AM rather than 
%%%    sortheader_importchan and allow for snipfile_overide argument
%%% -allow for dirname to be an absolute path if it contains current directory 
% Syntax:
%   cass_sort_apply(dirname,options)
% where
%   dirname is the name of a directory which contains the standard
%     structure, i.e., directories of the name chanX, where X is a
%     channel number; these chanX directories contain *.mat files with
%     sorting info.  Note that "autosort_info.mat," or any other file
%     which does not have a unique set of sorting choices, is omitted
%     from application.
%   options is a stucture which may have the following fields:
%     usecurrentpath (default true): if true, will ignore the information
%       saved in fh.abspathstr and instead replace it with the current
%       directory path (if 0 or absent, will rely on fh.abspathstr)
%     chanfile: if supplied (see PREPARE_CHANFILE for syntax), then
%       cass_sort_apply will be run only on the supplied list of channels &
%       sorting files.
%
% This creates new files with the same basename, and extension
% ".sorted". They are, in fact, .mat format files.  These files contain
% the sorting decisions made by LANDMARK_SORT_APPLY.
% 
% See also: LANDMARK_SORT_APPLY, EPHYS.
%
% HISTORY: 
%   2005-08-15 load statement altered to include "pwd" so that 
%              ghost saves can't feed incorrect information (RCH)
%   2005-11-07 added "options" option, with only current option being
%              options.usecurrentpath
%   2005-12-05 added ./ to definition of dirchan so as not to confuse
%              crappy FAT32 partition reader 
%   2006-06-12 switched default setting of "usecurrentpath" options to 1
  
% Copyright 2005 by Timothy E. Holy
  
  if nargin < 2
    options = struct;
  end
  options = default(options,'usecurrentpath',true);
  
  if isfield(options,'chanfile')
    channels = [options.chanfile.channel];
  else
    channels = sorted_get_channels(dirname);
  end
  if isempty(channels)
    error([dirname ' did not contain any channels']);
  end
  
  %% AM added 8/10/2015
  pwd_length = length(pwd);
  if length(dirname)>=pwd_length && strcmp(dirname(1:pwd_length),pwd)
      abspath = [];
  else
      abspath = [pwd filesep];
  end
  %%
  
  tfilecontents = load([abspath dirname filesep 'overview']); %% AM edited 8/10/2015
  
  if ~isfield(tfilecontents,'sorthead')
    error('Don''t recognize format of overview file')
  end
  sorthead = tfilecontents.sorthead;
  
  for i = 1:length(channels)
    currentChannel = channels(i);
    
    %% AM added 8/10/15
    if length(dirname)>=pwd_length && strcmp(dirname(1:pwd_length),pwd)
        parentpath = [];
    else
        parentpath = './';
    end
    %%
   
    dirchan = [parentpath dirname filesep 'chan' num2str(currentChannel) filesep];%% AM edited 8/10/2015
    
    if isfield(options,'chanfile')
      fls = {options.chanfile(i).file};
    else
      fls = dirbyname([dirchan '*.mat']);
      fls = setdiff(fls,'autosort_info.mat');
    end
    for j = 1:length(fls)
      tfilecontents = load([dirchan fls{j}]);
      if ~isfield(tfilecontents,'sort_info')
        continue
      end
      sort_info = tfilecontents.sort_info;
      if (length(sort_info) > 1 || length(sort_info.Rfactor) ~= 1)
        warning(['Skipping ' dirchan fls{j}])
        continue
      end
      % OK, this looks like a genuine file
      % Apply the landmark sorting
      if options.usecurrentpath == 1
        nDir = length(sorthead);
        for iDir = 1:nDir
          sorthead(iDir).fh.abspathstr = [cd filesep];
        end
      end
      shc = sortheader_importchan_AM(sorthead,currentChannel,snipfile_override);
      landmarkIndex = landmark_sort_apply(sort_info,shc);
      nfiles = length(shc);
      chanclust = cell(1,nfiles);
      for fileIndex = 1:nfiles
        clust = sort_info.landmarkClust(landmarkIndex{fileIndex});
        clabel = agglabel(clust+2);  % +2 for delete=-1
        clabel(1:min(2,length(clabel))) = []; % delete the delete & noise types
        for cellIndex = 1:length(clabel)
          chanclust{cellIndex,fileIndex} = ...
              double(shc(fileIndex).sniptimes(clabel{cellIndex}));
        end
      end
      [pathname,basename,ext] = fileparts([dirchan fls{j}]);
      outputfile = [basename '.sorted'];
      save([dirchan outputfile],'chanclust')
      fprintf('  Completed %s%s\n',dirchan,fls{j});
    end
  end
