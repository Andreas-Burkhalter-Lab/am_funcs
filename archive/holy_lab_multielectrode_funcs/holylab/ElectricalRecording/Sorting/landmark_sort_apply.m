function landmarkIndex = landmark_sort_apply(sort_info,shc)
% LANDMARK_SORT_APPLY: apply a landmark sort to all snippets
% Syntax:
%   landmarkIndex = landmark_sort_apply(sort_info,shc)
% where
%   sort_info is sorting data of the type described in save_sort_info;
%   shc is a sortheader structure array specialized to a particular
%     channel (with SORTHEADER_IMPORTCHAN);
% and
%   landmarkIndex is a cell array of vectors (one for each file in shc)
%     containing the index of the dclosest landmark to each spike snippet.
%
% See also: CASS_SORT_APPLY, AUTOSORT, CASS, SORTHEADER_IMPORTCHAN.
  
% Copyright 2005 by Timothy E. Holy

  nsnips_in_mem = 5000;
  nsnips_per_file = [shc.numofsnips];
  nfiles = length(shc);
  landmarkIndex = cell(1,nfiles);
  % Project the landmark waveforms, if necessary
  landmarkPosition = sort_info.landmarkWaveform;
  if sort_info.use_projection
    pd = sort_info.projectDirections';
    landmarkPosition = pd*landmarkPosition;
  end
  % Find out whether a proj.mat file is already present
  have_projections = false;
  if isfield(shc,'dirname')
    fullpath = [shc(i).dirname filesep 'chan' num2str(channel)];
    projfilename = fullfile(fullpath,'proj.mat');
    if exist(projfilename,'file')
      tt = load(projfilename,'snipProj');
      snipProj = tt.snipProj;
      have_projections = true;
    else
      error('I expected to find the proj.mat file, but it wasn''t there');
    end
  end
  % Append time to the landmarks
  if sort_info.t2V
    landmarkPosition = [landmarkPosition; sort_info.landmarkT*sort_info.t2V];
  end
  tstart = sortheader_absolute_starttime(shc); % File start time
  for idx = 1:nfiles
    landmarkIndex{idx} = nan(1,nsnips_per_file(idx));
    nblocks = ceil(nsnips_per_file(idx)/nsnips_in_mem);
    % Load snippets in blocks small enough to fit in memory
    for j = 1:nblocks
      startpos = (j-1)*nsnips_in_mem + 1;
      endpos = min(startpos+nsnips_in_mem-1,nsnips_per_file(idx));
      if have_projections
        snips = snipProj{idx}(:,startpos:endpos);
      else
        snips = sortheader_readsnips(shc(idx),startpos:endpos);
        if sort_info.use_projection
          % Project waveform
          snips = pd*snips;
        end
      end
      % Append time to the snippets
      if sort_info.t2V
        sniptimes = double(shc(idx).sniptimes(startpos:endpos)');
        absoluteSnipTime = sniptimes/shc(idx).scanrate ...
          + tstart(idx) - min(tstart);
        snips = [snips; absoluteSnipTime * sort_info.t2V];
      end
      % Find the closest landmark to each snippet
      [dist,landmarkIndex{idx}(startpos:endpos)] = ...
          mindist(snips,landmarkPosition);
    end
  end
  
  