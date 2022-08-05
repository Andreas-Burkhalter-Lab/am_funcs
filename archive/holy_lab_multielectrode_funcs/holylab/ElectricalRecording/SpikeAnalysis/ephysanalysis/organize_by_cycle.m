function cycle_index = organize_by_cycle(ephysin)
% ORGANIZE_BY_CYCLE: organize a set of stimulus presentations into cycles
% A cycle is defined as a set of trials which are "near-neighbors" in
% terms of time of occurrence, which test each of several stimuli.  The
% cycles do not need to test the stimuli in fixed order; for example,
% each cycle might present the fixed stimuli in random order.  However,
% the central notion is that each stimulus will be tested at most once
% per cycle, so this notion is not well-suited to designs in which
% stimuli are chosen completely randomly.
%
% Syntax:
%   cycle_index = organize_by_cycle(ephysin)
% where
%   ephysin is a (tagged) structure array of ephys structures, supplied
%     in temporal order of occurrence;
% and
%   cycle_index is a cell array, each element containing the indices of
%     elements of the ephys structure array comprising the given cycle.

% Copyright 2005 by Timothy E. Holy
  
  % Organize the trials by their tags
  [utags,tmp,utag_index] = unique({ephysin.tag});
  n_tags = length(utags);
  [clabel,nlabel] = agglabel(utag_index);
  % Compute the maximum number of repeats
  n_cycles = max(nlabel);
  % Determine where to put the breaks between cycles, in terms of the
  % earliest position needed to insure that two repeats of the same
  % stimulus do not end up in the same trial.
  break_index = inf(1,n_cycles);
  for i = 1:n_tags
    break_index(1:nlabel(i)) = min(clabel{i},break_index(1:nlabel(i)));
  end
  break_index = [break_index length(ephysin)+1];  % to simplify later stuff
  % Organize the trials by cycles
  cycle_index = cell(1,n_cycles);
  for i = 1:n_cycles
    cycle_index{i} = break_index(i):break_index(i+1)-1;
  end
  