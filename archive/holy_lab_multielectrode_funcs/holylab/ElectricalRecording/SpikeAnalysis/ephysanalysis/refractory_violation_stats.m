function outdata = refractory_violation_stats(ephys,params)
% REFRACTORY_VIOLATION_STATS: statistics on refractory violations
% Syntax:
%    rdata = refractory_violation_stats(ephys)
%    rdata = refractory_violation_stats(ephys,params)
% where
%   ephys is an ephys structure array with celltimes already loaded;
%   params is a structure that can have the following fields
%     refract_period_in_ms (default 20): The period below which a pair of
%       spikes is considered a refractory violation, supplied in
%       milliseconds
%     min_period_in_ms (default 0): You can exclude double-triggers and the
%       like by supplying a value larger than 0 (but smaller than the
%       refractory period). In milliseconds.
% and
%   rdata is a structure array, one element per cell, with the following
%   data:
%     cellnum: the "name" of the cell
%     n_spikes: the total number of spikes fired by this cell
%     n_violating_spikes: the total number of spikes violating the
%       refractory period
%     tac: the interspike intervals for all violations.

% Copyright 2007 by Timothy E. Holy

  if (nargin < 2)
    params = struct;
  end
  params = default(params,'refract_period_in_ms',20);
  params = default(params,'min_period_in_ms',0);
  
  n_intervals = length(ephys);
  ms2scans = [ephys.scanrate] / 1000;
  params.refract_period_in_scans = params.refract_period_in_ms * ms2scans;
  params.min_period_in_scans = params.min_period_in_ms * ms2scans;
    
  % Check to see that all elements contain info about the same cells
  cellnums = ephys(1).cellnums;
  for intervalIndex = 1:n_intervals
    if ~isequal(ephys(intervalIndex).cellnums,cellnums)
      error('The cell numbers are not consistent');
    end
  end
  
  n_cells = length(cellnums);
  for cellIndex = 1:n_cells
    tac = cell(1,n_intervals);
    n_spikes = zeros(1,n_intervals);
    for intervalIndex = 1:n_intervals
      trefr = params.refract_period_in_scans(intervalIndex);
      tmin = params.min_period_in_scans(intervalIndex);
      t = ephys(intervalIndex).celltimes{cellIndex};
      tac_tmp = autocorrspike(t,trefr);
      tac_tmp = tac_tmp(tac_tmp > tmin);
      tac{intervalIndex} = tac_tmp;
      n_spikes(intervalIndex) = length(t);
    end
    tac = cat(2,tac{:});
    outdata(cellIndex) = struct('cellnum',cellnums(cellIndex),...
				'n_spikes',sum(n_spikes),...
				'n_violating_spikes',length(tac),...
        'tac',tac);
  end