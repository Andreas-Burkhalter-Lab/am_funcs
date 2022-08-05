function inten = get_trials_dfof(intensities, onset_list, trange)
% GET_TRIALS_DFOF: a simple function to get temporal traces from imaging data
%
% Syntax:
%   inten = get_trials_dfof(intensities, onset_list, trange)
% where
%   intensities is a matrix or vector of intensity as a function of time
%     (if supplied as a matrix, is should be of size n_cells-by-n_stacks)
%   onset_list is a vector of valve onset times, or a cell array of onset
%     times for different values (e.g., see find_stimulus_start)
%   trange is a 2-vector, e.g., [-5 10], specifying the range of time
%     indices (inclusive) relative to valve onset.
% and
%   inten is a matrix of size n_cells-by-T-by-n_trials, where n_trials is
%     the # of trials and T is the duration of the timerange
%     (diff(trange)+1), containing the response on each trial. Each is
%     expressed as a deltaF/F, where the baseline period ranges from the
%     beginning of the trial's timerange to just before valve onset.
%   (if onset_list is a cell array, inten is instead a cell array of
%    matrices of the type described above)
%
% See also: GET_TRIALS_DFOF_CORRECTED, FIND_STIMULUS_START.

% Copyright 2010 by Timothy E. Holy & Diwakar Turaga

  [n_cells,n_stacks] = size(intensities);
  if (n_stacks == 1)
    intensities = intensities';
    n_cells = size(intensities,1);
  end
  output_matrix = false;
  if ~iscell(onset_list)
    onset_list = {onset_list};
    output_matrix = true;
  end
  inten = cell(1,length(onset_list));
  for i = 1:length(onset_list)
    inten{i} = zeros(n_cells,diff(trange)+1,length(onset_list{i}));
    for trial_indx = 1:length(onset_list{i})
      onset = onset_list{i}(trial_indx);
      inten_snip = intensities(:,onset+trange(1):onset+trange(2));
      base_line = mean(inten_snip(:,1:-trange(1)),2);
      inten_snip_norm = bsxfun(@rdivide,inten_snip,base_line) - 1;
      inten{i}(:,:,trial_indx) = inten_snip_norm;
    end
  end
  if output_matrix
    inten = inten{1};
  end
