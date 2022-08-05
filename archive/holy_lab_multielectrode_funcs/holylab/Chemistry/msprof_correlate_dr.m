function [lsqerr,clabel,mapc2dr] = msprof_correlate_dr(concdata,drdata,options)
% MSPROF_CORRELATE_DR: correlate physiological response with concentration
%
% Syntax:
%   [lsqerr,clabel,mapc2dr] = msprof_correlate_dr(concdata,drdata,options)
% where
%   concdata is a structure of the form output by msprof_concentrations;
%   drdata is a structure of the form described in msprof_monotonic_dr;
%   options is a structure which may have the following fields:
%     mode: 'dr_monotonic' (default) fits the firing rates to a model that
%       is monotonic with concentration; 'sensitivity_monotonic' fits the
%       sensitivity (the summed deltar across log-spaced concentrations) to
%       a monotonic model; 'hill' fits the responses to a Hill model. For
%       the latter, see fit_hill_equation for additional options.
%     sigma0ratio (default 1): the error bars can be mis-estimated
%       if there are many stimuli and few trials. This option allows you to
%       add a fraction of the rms error across stimuli to each error
%       estimate,
%          newerr = sqrt(olderr.^2 + options.sigma0ratio*mean(olderr.^2))
% and
%   lsqerr is an n_cells-by-n_compounds matrix, holding the fitting error;
%   clabel is a cell array of stimulus names;
%   mapc2dr is an index vector, matching the stimuli in drdata to
%     the stimuli in concdata. In particular,
%       drdata.stimlabel{mapc2dr(i)} = concdata.stimtag{i}
%     if mapc2dr(i) is not NaN.
%
% See also: MSPROF_CORRELATE_GUI, MSPROF_MONOTONIC_DR (which gives more detailed output).

% Copyright 2009-2010 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'sigma0ratio',1,'mode','dr_monotonic');
  
  % Use a shrinkage-like operation to regularize the errors: add in a
  % component of the mean error across trials
  rateerrs = drdata.rateerrs;
  mse = mean(rateerrs.^2,1);
  rateerrs = sqrt(rateerrs.^2 + options.sigma0ratio*repmat(mse,size(rateerrs,1),1));
  
  n_compounds = length(concdata.mz);
  % We have to match entries in the concdata to entries in drdata. Because
  % there are different numbers of each (due to multiple concentrations in
  % the drdata) this gets a bit tricky.
  % First find the common stimuli
  clabel = intersect(concdata.stimtag,drdata.uniquelabels);
  if isempty(clabel)
    error('There are no stimuli in common.')
  end
  % Now find which are matched to the common list
  flagdr = ismember(drdata.stimlabel,clabel);
  flagc = ismember(concdata.stimtag,clabel);
  % For each (matching) entry in drdata, find the corresponding entry in
  % concdata
  mapc2dr = nan(1,length(drdata.stimlabel));
  mapc2dr(flagdr) = findainb(drdata.stimlabel(flagdr),concdata.stimtag(flagc));
  % Now map the unique labels in drdata to concdata (needed for sensitivity
  % correlation)
  flagsens = ismember(drdata.uniquelabels,clabel);
  mapc2sens = nan(1,length(drdata.uniquelabels));
  mapc2sens(flagsens) = findainb(drdata.uniquelabels(flagsens),concdata.stimtag(flagc));
  % Extract matching elements of the firing rate
  rates = drdata.rates(flagdr,:);
  rateerrs = rateerrs(flagdr,:);
  concfac = drdata.stimconc(flagdr);
  n_cells = size(drdata.rates,2);
  lsqerr = zeros(n_cells,n_compounds);
  fprintf('%d cells to process: ',n_cells);
  for cellIndex = 1:n_cells
    fprintf('%d...',cellIndex);
    tmpr = rates(:,cellIndex);
    tmpre = rateerrs(:,cellIndex);
    for compoundIndex = 1:n_compounds
      tmpI0 = concdata.I(mapc2dr(flagdr),compoundIndex);
      tmpI = tmpI0 .* concfac(:); % Contains estimate of ligand concentration in each presented sample
      switch options.mode
        case 'dr_monotonic'
          [~,sortOrder] = sort(tmpI);  % Sort them into increasing concentration
          [~,lsqerr(cellIndex,compoundIndex)] = lsq_monotonic(tmpr(sortOrder),tmpre(sortOrder)); % Solve for the idealized firing rate and total mismatch
        case 'sensitivity_monotonic'
          relconc = drdata.relconc(flagsens,:);
          relconcerr = drdata.relconcerr(flagsens,:);
          tmpI = concdata.I(mapc2sens(flagsens),compoundIndex);
          [~,sortOrder] = sort(tmpI);  % Sort them into increasing concentration
          [~,lsqerr(cellIndex,compoundIndex)] = lsq_monotonic(relconc(sortOrder,cellIndex),relconcerr(sortOrder,cellIndex));
        case 'hill'
          res = fit_hill_equation(tmpI,tmpr,tmpre,options);
          if ~isnan(res.full.chisq)
            lsqerr(cellIndex,compoundIndex) = res.full.chisq;
          elseif ~isnan(res.linear.chisq)
            lsqerr(cellIndex,compoundIndex) = res.linear.chisq;
          else
            lsqerr(cellIndex,compoundIndex) = res.flat.chisq;
          end
        otherwise
          error('Mode %s not recognized',options.mode);
      end
    end
  end
  fprintf('\n');

  