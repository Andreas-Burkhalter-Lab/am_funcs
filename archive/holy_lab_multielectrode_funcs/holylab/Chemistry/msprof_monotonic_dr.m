function [C,R,err,stimOrder] = msprof_monotonic_dr(concdata,drdata)
% MSPROF_MONOTONIC_DR: compare concentration & firing rate data using a monotonic model
%
% This function takes a set of MS profile data across samples, and compares
% physiological responses to the same samples (which may be presented at a
% variety of concentrations). For each compound, it computes the
% concentration in each sample presented physiologically. It then fits a
% firing rate curve whose only constraint is that it be monotonic in
% concentration.
%   
% Syntax:
%  [C,R,err,stimOrder] = msprof_monotonic_dr(concdata,drdata)
% where
%   concdata has fields:
%     stimlabel: a 1-by-n_samples cell array containing a string label for
%       each type of stimulus
%     I: an n_samples-by-n_compounds matrix containing the total ion count
%       for each compound
%   drdata has fields:
%     stimlabel: a 1-by-n_stimuli cell array containing string labels for
%       each type of stimulus (note concentration information should _not_
%       be part of the string, just identity). These should match those in
%       concdata.stimlabel, but they do not have to be in the same order or
%       even contain all of the same entries.
%     stimconc: a 1-by-n_stimuli vector containing the relative
%       concentration at which each stimulus was presented
%     rates and rateerrs: n_stimuli-by-n_cells matrices, giving the mean
%       and s.e.m. across trials.
% and
%   C is an n_compounds-by-n_stimuli matrix containing each stimulus'
%     concentration of the given compound
%   R is an n_compounds-by-n_stimuli-by-n_cells array containing the
%     best-fit firing rates from the monotonic model
%   err is an n_compounds-by-n_cells matrix containing the fitting error of
%     the monotonic model
%   stimOrder is a vector indexing the stimuli in drdata in the order of C
%     and R.
%
% See also: MSPROF_CORRELATE_DR.

% Copyright 2009 by Timothy E. Holy

  % Match the stimuli between concdata and drdata
  [ul,tmp,stimIndex] = unique(drdata.stimlabel);
  [stimGroup,nstimGroup] = agglabel(stimIndex);
  [stimlabel,indexc,indexdr] = intersect(concdata.stimlabel,ul);
  stimOrder = cat(2,stimGroup{indexdr});
  stimIndex = stimIndex(stimOrder);
  stimConc = drdata.stimconc(stimOrder);
  stimConc = stimConc(:);
  % Prepare the output
  n_compounds = size(concdata.I,2);
  n_stim = sum(nstimGroup(indexdr));
  n_cells = size(drdata.relconc,2);
  R = zeros(n_compounds,n_stim,n_cells);
  C = zeros(n_compounds,n_stim);
  err = zeros(n_compounds,n_cells);
  for compoundIndex = 1:n_compounds
    % Compute the concentration of this species in each sample
    I = concdata.I(indexc,compoundIndex);
    c = I(stimIndex) .* stimConc;
    C(compoundIndex,:) = c';
    % Sort the samples in order of increasing concentration
    [cs,sortOrder] = sort(c);
    % Collect the firing rate data
    for cellIndex = 1:n_cells
      r = drdata.rates(stimOrder,cellIndex);
      sigma = drdata.rateerrs(stimOrder,cellIndex);
      % Fit the firing rate data with a model that is monotonic in
      % concentration
      [Rtmp,err(compoundIndex,cellIndex)] = ...
        lsq_monotonic(r(sortOrder),sigma(sortOrder));
      R(compoundIndex,sortOrder,cellIndex) = Rtmp;
    end
  end
  % Prepare the last elements of the output
  stimconc = stimConc';
  stimlabel = ul(stimIndex);
  