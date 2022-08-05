function [t2V,options] = autosort_calc_t2V(shc,options)
% AUTOSORT_CALC_T2V: calculate the tradeoff between time and voltage
%
% This function tries to assess how quickly spike waveforms are changing
% over time, and measures that rate in terms of a constant which converts
% times to voltages.  It performs this task in the following way:
%   1.  It loads several distinct blocks of spikes, separated in time.
%   2.  It creates a set of landmarks (using K-means), and assigns all
%       spikes to a landmark.
%   3.  It measures the consistency of this assignment between adjacent
%       blocks, using the algorithm in ESTIMATE_CLUST_MOVEMENT to
%       estimate the typical length scale of changes between blocks.
%   4.  It converts these lengths to a rate by dividing by the interval
%       between blocks, and then takes the median value of ell/t across
%       adjacent block pairs.
%
% Syntax:
%   t2V = autosort_calc_t2V(shc,options)
%   [t2V,options] = autosort_calc_t2V(shc,options)
% where
%   shc is a sortheader specialized to a specific channel, see
%     SORTHEAD_IMPORTCHAN;
%   options is a structure which may have the following fields:
%      t2V_n_blocks: the number of blocks of spikes used in measuring the
%        time rate of change of waveform cluster positions. Default: 4
%      t2V_min_n_blocks: the minimum number of blocks considered
%        acceptable for calculating t2V. Default: 2
%      t2V_maxblocksize: the maximum number of spike waveforms in each
%        block.  Default: 1000
%      t2V_minblocksize: the minimum number of spike waveforms in each
%        block.  Default: 200
%      t2V_n_landmarks: the number of landmarks to use.  Default: 20
% and
%   t2V is a scalar, the sought-after conversion factor (t in units of
%     seconds);
%   options is an output options structure with any unsupplied fields
%     filled in with their default values.
%
% See also: ESTIMATE_CLUST_MOVEMENT, SORTHEAD_IMPORTCHAN.
 
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  options = autosort_calc_t2V_options(options);
  
  % Load several "chunks" of spikes, and then measure the
  % (in)consistency in pseudo-cluster positions between blocks
  % (pseudo-clusters = "clusters" identified by kmeans)
  [fileIndex,snipIndex] = choose_chunks(shc,options);
  if isempty(fileIndex)
    % Didn't have enough snippets for this analysis, abort
    t2V = 0;
    return
  end
  snips = {};
  tblock = [];
  nperblock = [];
  for blockIndex = 1:length(fileIndex)
    % Load block of snippets
    snips{blockIndex} = sortheader_readsnips(shc(fileIndex(blockIndex)),...
                                             snipIndex{blockIndex});
    % Determine the time of the first snippet in each block
    tblock(blockIndex) = ...
        sortheader_absolute_sniptime(shc(fileIndex(blockIndex)),...
                                     snipIndex{blockIndex}(1));
    nperblock(blockIndex) = length(snipIndex{blockIndex});
  end
  snipsall = cat(2,snips{:});
  % Pseudo-cluster the snippets, and then compare pseudo-cluster
  % occupancy across blocks
  [pointLabel,center,rms_dist] = kmeans_hard(snipsall,options.t2V_n_landmarks);
  ell = estimate_clust_movement(pointLabel,rms_dist,nperblock);
  dt = diff(tblock);
  % And, put it all together
  t2V = median(ell/dt);

  
function options = autosort_calc_t2V_options(options)
  if ~isfield(options,'t2V_n_blocks')
    options.t2V_n_blocks = 4;
  end
  if ~isfield(options,'t2V_min_n_blocks')
    options.t2V_min_n_blocks = 2;
  end
  if ~isfield(options,'t2V_maxblocksize')
    options.t2V_maxblocksize = 1000;
  end
  if ~isfield(options,'t2V_minblocksize')
    options.t2V_minblocksize = 200;
  end
  if ~isfield(options,'t2V_n_landmarks')
    options.t2V_n_landmarks = 20;
  end

function [fileIndex,snipIndex] = choose_chunks(shc,options)
% Find "chunks" of spikes within recordings; try to distribute these
% chunks as far apart as possible (definitely non-overlapping)
% This is surprisingly tricky!

  nFiles = length(shc);
  nSnipsPerFile = [shc.numofsnips];
  % Can we lay down even the minimum-sized blocks?
  curSize = options.t2V_minblocksize;
  nChunksPerFile = floor(nSnipsPerFile/curSize);
  if (sum(nChunksPerFile) < options.t2V_min_n_blocks)
    % Can't do it
    fileIndex = [];
    snipIndex = {};
    return;
  end
  % Do we have too many?
  newSize = curSize;
  while (sum(nChunksPerFile) > options.t2V_n_blocks && ...
         newSize < options.t2V_maxblocksize)
    % Grow the blocks
    curSize = newSize;
    newSize = min(round(1.5*curSize),options.t2V_maxblocksize);
    nChunksPerFile = floor(nSnipsPerFile/newSize);
  end
  % If the last one tried worked, keep it
  if (sum(nChunksPerFile) >= options.t2V_min_n_blocks && ...
      newSize == options.t2V_maxblocksize)
    curSize = newSize;
  end
  nChunksPerFile = floor(nSnipsPerFile/curSize);
  startIndex = cell(1,nFiles);
  fileIndex = cell(1,nFiles);
  for i = 1:nFiles
    startIndex{i} = ((1:nChunksPerFile(i))-1)*curSize+1;
    fileIndex{i} = i+zeros(1,nChunksPerFile(i));
  end
  startIndex = cat(2,startIndex{:});
  fileIndex = cat(2,fileIndex{:});
  n_blocks = length(startIndex);
  % If we still have too many, we have to pitch some
  if (n_blocks > options.t2V_n_blocks)
    keepIndex = linspace(1,n_blocks,options.t2V_n_blocks);
    keepIndex = round(keepIndex);
    startIndex = startIndex(keepIndex);
    fileIndex = fileIndex(keepIndex);
  end
  % Set up the final choice
  snipIndex = cell(1,length(fileIndex));
  for i = 1:length(snipIndex)
    snipIndex{i} = startIndex(i):startIndex(i)+curSize-1;
  end
