function save_sort_info(sort_info,filename)
% SAVE_SORT_INFO: save spike sorting information
%
% Note this is not actually used to save decisions about how individual
% spikes are assigned; instead, this saves the information needed to make
% such decisions.
%
% Syntax:
%   save_sort_info(sort_info, filename)
% where channel_sort_info is a structure with the following fields:
%   channel: an integer channel number;
%   sniprange: a two vector, [snipbegin snipend] giving the offset
%     relative to the peak of the spike (e.g., [-10 30]);
%   use_projection: a flag, if true then the waveforms are to be
%     projected to reduce dimensionality;
%   projectDirections: the directions (filters) used in projecting the
%     waveforms (present only if use_projection is true), each direction
%     being a column of this matrix;
%   t2V: the constant of proportionality converting times to voltages,
%     so that time can be used as an additional coordinate in sorting
%     (=0 if not using time as a coordinate);
%   landmarkWaveform: the spike waveform associated with each landmark
%     (in original voltage units, not projected by filters). This will be
%     a matrix, organized in column format (each column associated with
%     one landmark).
%   landmarkT: the time coordinate associated with each landmark, in
%     whatever units t2V is expressed in (will be a vector of length
%     n_landmarks);
% The next two fields are optional; if you don't supply them, you are
% basically providing the info needed for clustering without storing any
% of the merge results.
%   Rfactor: a scalar or vector of smoothing factors that were used (use
%     a scalar NaN if you want to encode that this is meaningless, e.g.,
%     for purely manual assignment);
%   landmarkClust: a vector of cluster numbers assigned to each landmark,
%     if Rfactors is a scalar; a matrix of cluster numbers, of size
%     n_Rfactors by n_landmarks, if Rfactors is a vector.
% The next is present only if some landmarks are assigned to the delete
% pile, and the user later tweaks the projectDirections:
%   projectDirDelete: the projectDirections to be used in deciding
%     whether a given waveform should be deleted or retained.  If a
%     snippet is retained, then its actual classification is determined
%     in terms of projectDirections. This field is here to acknowledge
%     the problem that spike artifacts can distort the projectDirections;
%     we might initially use one set of directions to identify artifacts
%     separately from spikes & noise, and another set of directions to
%     discriminate among spikes/noise.
%
% Note that sort_info can be a structure array, if you want to save
% replicates to a file.
%
% filename is an optional filename; if not supplied, it will first seek a
% subdirectory with the name "chanX", where X is a string corresponding
% to the channel number.  It will then save a file of name
% "sort_info.mat", in that subdirectory if found, or in the current
% directory otherwise.
%
% See also: LOAD_SORT_INFO, AUTOSORT, CASS.

% Copyright 2005 by Timothy E. Holy

  channel = unique([sort_info.channel]);
  chanStr = ['chan' num2str(channel)];
  [savedir,this_file] = fileparts(filename);
  if exist(savedir,'dir')
    savedir = '';   % we'll use the full path-qualified version passed in
  elseif exist([pwd filesep chanStr],'dir')
    savedir = [chanStr filesep];
  else
    savedir = '';
  end
  if (nargin < 2)
    filename = 'sort_info';
  end
  % Check it for valid syntax
  required_fieldnames = {'channel','sniprange','use_projection',...
                      't2V','landmarkWaveform','landmarkT'};
  missing_fieldnames = setdiff(required_fieldnames,fieldnames(sort_info));
  if ~isempty(missing_fieldnames)
    error('The required fieldname %s is not present in sort_info\n',...
          missing_fieldnames{:});
  end
  save([savedir filename],'sort_info')
