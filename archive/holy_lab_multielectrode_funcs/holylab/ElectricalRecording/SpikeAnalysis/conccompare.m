function results = conccompare(rates,errs,stimlabel,stimconc,options)
% CONCCOMPARE: estimate relative concentrations of stimulus components from deltar
% Given firing rates of a population of neurons to a complex stimulus
% blend, attempt to estimate the relative concentration of the particular
% ligand driving a particular neuron across the different stimuli.  At
% least one stimulus must be presented at two or more different
% concentrations, and more typically all stimuli will be presented at two
% or more concentrations.
%
% This function assumes classic first-order steady-state binding affinity
% of a single ligand to a single receptor, but of course may be used as
% an approximate measure for multiunit activity or in cases where more
% than one ligand might activate a given cell or recording site.  This
% function returns a goodness of fit so that the accuracy of the model
% can be evaluated.
%
% The fitting process goes something like this: for a given cell, find
% the "best stimulus" for this cell.  For this stimulus/cell pair, a real
% association constant will be fit.  The effectiveness of all other
% stimuli will be expressed in terms of a concentration relative to this
% "best stimulus."  This arrangement proves to be a good choice for
% factoring out some of the degeneracies in the parameters.
% 
% Syntax:
%   results = conccompare(rates,errs,stimlabel,stimconc,options)
%   results = conccompare(resultsold,options)
% where
%   rates is an nstimuli-by-ncells matrix of average steady-state firing
%     rate elevations (or deltars in a crude case);
%   errs is an nstimuli-by-ncells matrix of standard errors in the
%     rates;
%   stimlabel is a cell array of nstimuli strings, giving the identity of
%     each stimulus (with nunique different strings);
%   stimconc is a vector of size 1-by-nstimuli, giving the concentrations
%     of the individual stimuli;
%   options is an optional structure with the following fields:
%     plot: if true, graphical output is produced (you can also set any
%       options from CONCCOMPAREPLOT as well);
%     variablebaseline: if true, adds a fitting parameter to account for a
%       nonzero deltar independent of stimulus identity & concentration
%       (default false);
%     noclean: if true, skip the "cleaning" step in which parameters are
%       tested to see if the fit is genuine or simply an upper bound;
%     responsivethresh: responsiveness of the cell as a whole is tested
%       by making sure that an increase in all of the association
%       constants gives a noticeably worse fit.  This parameter
%       determines the number of chisq units by which a "real" response
%       should deteriorate (default: 3);
%     boundthresh: responsiveness to a particular stimulus is judged by
%       the same shifting mechanism as above for "responsivethresh,"
%       except that the association constant/relative concentration for
%       each stimulus is varied independently.  If a cell is judged as
%       not being responsive to a particular stimulus, an effective upper
%       bound for the stimulus is computed. boundthresh is also in chisq
%       units (default: 2);
%     shiftthresh: if a given stimulus does not drive a cell (by the
%       criterion above for "boundthresh"), this parameter determines the
%       increase in chisq which can be tolerated in shifting the curve to
%       the right (default: 1);
%     maxeval (default 500): maximum number of function evaluations in
%       minimization;
% and
%   results is structure, with the following fields:
%     options: a copy of the input options;
%     uniquelabels: a cell array of unique stimulus names;
%     rates, rateerrs, stimlabel, stimconc: copies of the input
%       information as described above;
%     dof: number of degrees of freedom in the fit (useful for
%       interpreting the significance of chisq);
%     beststim: a 1-by-ncells cell array of stimulus names, the best one
%       for each cell;
%     fitlabels: a nparameters-by-ncells cell array of names of the
%       different fitting parameters (references to particular stimuli
%       refer to their relative concentration);
%     fitparams: the best-fit parameters (prior to any
%       tweaking/validation);
%     chisq: the chisq for the fit;
%     invcov: inverses of the covariance matrix;
%     minimumfound: a 1-by-ncells vector indicating whether a fitting
%       minimum was found, or whether the fitting routine had to give up;
%     responsive: a 1-by-ncells vector, with 1 indicating a cell which
%       responds to some stimulus (not available if "noclean" is
%       specified);
%     responsivetest: the increment in chisq during the responsive
%       test. A cell is deemed responsive only if this increment is less
%       than the threshold determined in the input options structure;
%     rmax is a 1-by-ncells vector of maximal firing rates for the given
%       cell;
%     logK is the logarithm of the association constant for the "best
%       stimulus";
%     cOverK is a 1-by-ncells vector containing the estimate of c/K, the
%       concentration divided by the association constant, for the most
%       potent stimulus for the given cell;
%     relconc is an nunique-by-ncells matrix of estimated relative
%       concentrations of the ligand, with the most potent stimulus for a
%       given cell set at 1;
%     bound is an nunique-by-ncells matrix, where 1 is used whenever a
%       particular relative concentration is merely an upper bound (not
%       available if "noclean" is specified);
%     boundtest is the increment in chisq arising from the upper-bound test.
%
% See also: fit_hill_equation, CONCCOMPAREPLOT, CONCRATIO.
  
  if (nargin < 4 & isstruct(rates))
    if (nargin == 2)
      options = errs;
    else
      options = struct;
    end
    stimconc = rates.stimconc;
    stimlabel = rates.stimlabel;
    errs = rates.rateerrs;
    rates = rates.rates;
  elseif (nargin < 5)
      options = struct;
  end
  options = ccoptions(options);
  [nstimuli,ncells] = size(rates);
  if any(size(rates) ~= size(errs))
    error(['The dimensions on the rates does not agree with the dimensions' ...
           ' on the errors']);
  end
  if (~iscell(stimlabel) | length(stimlabel) ~= nstimuli)
    error('The stimulus labels do not agree with the number of stimuli');
  end
  % Create the fields of results, so they will always be in the same order
  % (that way can build an array across experiments)
  fn = {'options','uniquelabels','rates','rateerrs','stimlabel','stimconc',...
      'dof','beststim','fitlabels','fitparams','chisq','invcov','minimumfound',...
      'responsive','bound','responsivetest','rmax','logK','cOverK','relconc',...
      'boundtest','fitfailed'};
  cellfields = [8 9];
  for i = 1:length(fn)
    results.(fn{i}) = [];
  end
  for i = 1:length(cellfields)
    results.(fn{cellfields(i)}) = {};
  end
  results.options = options;  % So fitting options are recorded for posterity
  % Find the unique stimuli and the highest concentration of each
  uniquestim = unique(stimlabel);
  nunique = length(uniquestim);
  % Copy useful data to output
  results.uniquelabels = uniquestim;
  results.rates = rates;
  results.rateerrs = errs;
  results.stimlabel = stimlabel;
  results.stimconc = stimconc;
  for i = 1:nunique
    indx = strmatch(uniquestim{i},stimlabel,'exact');
    [tmp,maxcindx(i)] = max(stimconc(indx));
    maxcindx(i) = indx(maxcindx(i));
  end
  % Now find the best stimulus for each neuron
  % This is just done by finding the stimulus which maximally stimulates
  % the neuron (because of error bars, this could prove to be the wrong
  % choice upon further analysis, but we can correct it later)
  [maxr,bestindxu] = max(rates(maxcindx,:)); % This is indexed into unique stimuli
  bestindx = maxcindx(bestindxu); % This is indexed into stimuli
  results.dof = nstimuli - (nunique+1);
  for i = 1:length(bestindxu)
      results.beststim(i) = uniquestim(bestindxu(i));
  end
  % Now loop over individual neurons
  for i = 1:ncells
    % Set up initial guesses and massage parameters into a form suitable
    % for optimization
    params.rates = rates(:,i);
    params.rateerrs = errs(:,i);
    if any(params.rateerrs == 0)
      warning('sigma is zero, fixing')
      params.rateerrs(params.rateerrs == 0) = median(params.rateerrs);
    end
    if any(params.rateerrs == 0)
      warning('couldn''t fix, continuing onto the next cell')
      results.fitfailed(i) = true;
      continue
    end
    params.stimlabel = stimlabel;
    params.stimconc = stimconc;
    params.Xlabel = uniquestim(bestindxu(i));
    params.Xlabel = [params.Xlabel,setdiff(uniquestim,params.Xlabel{1})];
    X = [];
    X(1) = maxr(i);
    %X(1) = 5;
    X(2) = 0;
    %X(2) = log(0.02);
    for j = 1:nunique-1
      indx = strmatch(params.Xlabel{j+1},uniquestim,'exact');
      % This guesses ratios of concentrations as being the ratios of
      % firing rates
      X(2+j) = log(max([1 params.rates(maxcindx(indx))])/X(1));
    end
    if options.variablebaseline
      X(end+1) = 0;
      params.variablebaseline = 1;
    end
    % Do the optimization
    results.fitfailed(i) = false;
    try
      [Xbest,chisq,exitflag] = fminunc(@concchisq,X,optimset('LargeScale','off','MaxFunEvals',options.maxeval,'Display',options.display),params);
    catch ME
      fprintf('Error on i = %d\n',i);
      ME
      ME.stack
      results.fitfailed(i) = true;
      continue
    end
    results.fitlabels(:,i) = [{'maxr','logKbest'} params.Xlabel(2:end)];
    results.fitparams(:,i) = Xbest';
    results.chisq(i) = chisq;
    % Find the covariance
    [chisq,gradchisq,invcov] = concchisq(Xbest,params);
    results.invcov(:,:,i) = invcov;
    % Probe to see whether particular relative concentrations are real,
    % meaning that they cannot be decreased with impunity.  If they are not
    % real, determine the upper bound and indicate as much.
    results.minimumfound(i) = exitflag;
    if (options.noclean == 0)
        % First do all stimuli together: increase K, thereby decreasing the
        % effective concentration of all stimuli. If this doesn't increase
        % chi-squared significantly, then there wasn't a real response to
        % any stimulus.
        results.responsive(i) = 1;
        results.bound(1:nunique,i) = zeros(nunique,1);
        chi2 = concchisq(Xbest + [0 log(10) zeros(1,length(Xbest)-2)],params);
        results.responsivetest(i) = chi2-chisq;
        if (chi2-chisq < options.responsivethresh*results.dof)
            % This whole cell is bogus
            results.responsive(i) = 0;
        else
            % Now do all but the best stimulus (the best is covered by the
            % previous step). Decrease the relative concentration of the
            % stimulus by tenfold; if this doesn't hurt chi-squared, then
            % this stimulus was not effective in exciting the cell.
            % In this case, boost its relative concentration to the maximum
            % allowed by the data, and then indicate that it's an upper
            % bound.
            Xbestnew = Xbest;
            for j = 1:nunique-1
                indx = strmatch(params.Xlabel{j+1},uniquestim,'exact');              
                Xmanip = [zeros(1,j+1) 1 zeros(1,length(Xbest)-j-2)];
                chi2 = concchisq(Xbest - log(10)*Xmanip,params);
                results.boundtest(j,i) = chi2-chisq;
                if (chi2-chisq < options.boundthresh)
                    % It's not a real fit; the best we can get is an upper
                    % bound. Increase by 10% each time.
                    results.bound(indx,i) = 1;
                    chi2 = concchisq(Xbest + log(1.1)*Xmanip,params);
                    Xbesttemp = Xbest;
                    while (chi2 - chisq < options.shiftthresh)
                        Xbesttemp = Xbesttemp + log(1.1)*Xmanip;
                        chi2 = concchisq(Xbesttemp + log(1.1)*Xmanip,params);
                    end
                    Xbestnew(j+2) = Xbesttemp(j+2);
                end
            end
            Xbest = Xbestnew;
        end
    end
    % Now massage parameters into sensible outputs
    results.rmax(i) = Xbest(1);
    results.logK(i) = Xbest(2);
    results.cOverK(i) = stimconc(bestindx(i))/exp(Xbest(2));
    relconc1 = exp(Xbest(3:nunique+1));
    results.relconc(:,i) = [relconc1(1:bestindxu(i)-1), 1, relconc1(bestindxu(i):end)]';
    if options.variablebaseline
      results.r0(i) = Xbest(end);
    end
    % Check to make sure that the maximum relative concentration is 1, and
    % make adjustments as needed
    %if any(results.relconc(:,i) > 1)
    %    rcmax = max(results.relconc(:,i));
    %    results.relconc(:,i) = results.relconc(:,i)/rcmax;
    %    results.cOverK(i) = results.cOverK(i)*rcmax;
    %end
  end
  % If desired, plot the results
  if (isfield(options,'plot') & options.plot)
      if ~isfield(options,'convergeonly')
          options.convergeonly = 1;
      end
      conccompareplot(results,options)
  end
  
  % Make sure check for relconc > 1, and adjust cOverK, relconc, etc.

  
  function options = ccoptions(options)
  if ~isfield(options,'noclean')
    options.noclean = 0;
  end
  if ~isfield(options,'plot')
    options.plot = 0;
  end
  if ~isfield(options,'responsivethresh')
    options.responsivethresh = 3;
  end
  if ~isfield(options,'boundthresh')
    options.boundthresh = 2;
  end
  if ~isfield(options,'shiftthresh')
    options.shiftthresh = 1;
  end
  if ~isfield(options,'variablebaseline')
    options.variablebaseline = 0;
  end
  if ~isfield(options,'maxeval')
    options.maxeval = 500;
  end
  if ~isfield(options,'display')
    options.display = 'final';
  end