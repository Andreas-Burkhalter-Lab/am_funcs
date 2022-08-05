function Icorr = intensity_corrected(I,onset_list,trange,b,Ictrl)
% intensity_corrected: subtract control responses from intensities vs. time
%
% This function allows you to generate a "Ringer's subtracted" full time
% course of intensities.
% 
% Syntax:
%   Icorr = intensity_corrected(I,onset_list,trange,b,Ictrl)
% where
%   I is a maxtrix or vector of intensity as a function of time (if
%     supplied as a matrix, is should be of size n_cells-by-n_stacks)
%   onset_list is a cell array, each element containing a vector of valve
%     onset times for each stimulus (see find_stimulus_start)
%   trange is a 2-vector, e.g., [-5 10], specifying the range of time
%     indices (inclusive) relative to the times in onset_list.
%   b is a n_stimuli-by-n_trials matrix, containing the amplitude of the
%     control response subtracted from each trial (in units of the mean
%     control response).
%   Ictrl is the response to the control stimulus.
% and
%   Icorr is the corrected intensity matrix.
%
% Example:
%  % (Assume you have "header" and "intensities" defined)
%  onset_list = find_stimulus_start(header.stim_lookup);
%  ctrlindex = 2;
%  trange = [-4 14];
%  [I,b,Ictrl] = get_trials_dfof_corrected(intensities,onset_list,trange,ctrlindex);
%  Icorr = intensity_corrected(intensities,onset_list,trange,b,Ictrl);

% Copyright 2010 by Timothy E. Holy

  Icorr = I;
  n_stimuli = length(onset_list);
  for stimIndex = 1:n_stimuli
    n_trials = length(onset_list{stimIndex});
    for trialIndex = 1:n_trials
      t = onset_list{stimIndex}(trialIndex);
      bg = mean(Icorr(:,t+trange(1):t-1),2);
      trng = t+trange(1):t+trange(2);
      % Compute the deltaF
      df = b{stimIndex}(trialIndex) * bsxfun(@times,bg,Ictrl);
      % Subtrace this deltaF from the intensity
      Icorr(:,trng) = Icorr(:,trng) - df;
    end
  end
  
      