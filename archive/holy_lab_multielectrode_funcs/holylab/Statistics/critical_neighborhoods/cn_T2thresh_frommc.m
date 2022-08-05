function [T2thresh,T2lsthresh] = cn_T2thresh_frommc(d,nmax,params)
  % cn_T2thresh_frommc: interpolate T^2 thresholds from Monte Carlo results
  % Syntax:
  %   [T2thresh,T2lsthresh] = cn_T2thresh_frommc(d,nmax,params)
  % where
  %   d is the desired dimensionality
  %   nmax is the total number of points in your data set
  %   params is a structure with fields 'covarianceModel' and 'pvalue'
  % and
  %   T2thresh is a vector containing the T^2 threshold for neighborhoods
  %     of sizes 1:nmax, from a fixed basepoint
  %   T2lsthresh is the same quantity for use during line search.
  %
  % This function relies on having Monte Carlo results stored to disk.
  %
  % See also:  cn_mclinesearch, cn_run_mc.
  
  % Copyright 2012 by Timothy E. Holy
  
  % Load the results from the Monte Carlo simulation
  switch params.covarianceModel
    case 'isotropic'
      load('T2iso');
    case 'diagonal'
      load('T2diag');
    case 'full'
      load('T2full');
    otherwise
      error('covarianceModel not recognized');
  end
  % Check bounds
  if (d > max(T2data.dlist))
    error('Cannot interpolate, dimensionality is larger than run in the simulation');
  end
  if (params.pvalue < min(T2data.pvaluelist) || params.pvalue > max(T2data.pvaluelist))
    error('Cannot interpolate, the p-value is outside the range of the simulation');
  end
  ninterp = 1:min(nmax,max(T2data.nlist));
  % Build the arrays that specify the lookup values. We'll use
  % interpolation on the logarithm of the parameters, so we take the log
  % first
  X = cell(1,3);
  [X{:}] = ndgrid(log(T2data.dlist),log(T2data.nlist),log(T2data.pvaluelist));
  % Do the interpolation
  T2thresh = interpn(X{:},T2data.basepoint,log(d),log(ninterp),log(params.pvalue));
  T2lsthresh = interpn(X{:},T2data.linesearch,log(d),log(ninterp),log(params.pvalue));
  % Pad with the last value up to nmax
  T2thresh(end+1:nmax) = T2thresh(end);
  T2lsthresh(end+1:nmax) = T2lsthresh(end);
end

  