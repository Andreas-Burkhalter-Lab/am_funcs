function [onset_list,ustim] = find_stimulus_start(stim_lookup,flushvalve,trange)
% FIND_STIMULUS_START: locate the first frame, organized by stimulus presentation
% Syntax:
%   [onset_list,ustim] = find_stimulus_start(stim_lookup)
%   [onset_list,ustim] = find_stimulus_start(stim_lookup,flushvalve)
%   [onset_list,ustim] = find_stimulus_start(stim_lookup,flushvalve,trange)
% where
%   stim_lookup is a vector containing a # associated with the stimulus
%     presented during each stack;
%   flushvalve is the valve number assigned to the flush (if missing or
%     empty: defaults to the most commonly-presented valve)
%   trange, if present, is a 2-vector of the form [first last]; any trials
%     that begin before first, or after last, are excluded;
% and
%   onset_list is a cell array, one element for each unique stimulus
%     number. Each element of this array contains the list of stack
%     numbers at which this stimulus was turned on.
%   ustim is a list of unique stimulus #s. The sequence of these
%     corresponds to the elements of onset_list.
%
% The second syntax, where you specify 'flushvalve', allows you to
% explicitly indicate the # corresponding to the flush valve. (By default
% it uses the most common valve # for the flush.)

% Copyright 2010 by Diwakar Turaga & Timothy E. Holy
  
  if (nargin < 2 || isempty(flushvalve))
    flushvalve = mode(stim_lookup);
  end
  if (nargin < 3 || isempty(trange))
    trange = [1 length(stim_lookup)];
  end
  % Find the stack # at which each stimulus turns on
  onsetIndex = find(stim_lookup(1:end-1) == flushvalve & ...
    stim_lookup(2:end) ~= flushvalve)+1;
  % Find the list of unique stimulus #s (excluding flush)
  [ustim,tmp,uIndex] = unique(stim_lookup(onsetIndex)); %#ok<ASGLU>
  % Collect the cases corresponding to the same stimulus
  onset_list = agglabel(uIndex);
  for i = 1:length(onset_list)
    tmp = onsetIndex(onset_list{i});
    tmp(tmp < trange(1)) = [];
    tmp(tmp > trange(2)) = [];
    onset_list{i} = tmp;
  end
  
