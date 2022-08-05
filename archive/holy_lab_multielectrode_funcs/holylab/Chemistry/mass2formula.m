function [formula,mass_err] = mass2formula(mass,isotopeinfo,options)
% MASS2FORMULA: generate list of possible molecular formulas producing a given mass
% Syntax:
%   formula = mass2formula(mass,isotopeinfo)
%   [formula,mass_error] = mass2formula(mass,isotopeinfo,options)
% where
%   mass is a scalar giving the desired mass (typically with 3-5 decimal
%     places of precision);
%   isotopeinfo is the output of ISOTOPES;
%   options may have the following fields:
%     mass_fractonal_error (default 1e-5): the amount of error, expressed
%       in terms of the fraction of the mass, that can be tolerated and still
%       considered as a match;
%     n_to_try_max (default Inf): use this to prevent long computation
%       times, will cause the routine to exit without even starting if the
%       number of possibilities is very large.  The number of different
%       elements is the primary determinant of the number of possible
%       combinations, so if needed call isotopes with fewer elements.
% and
%   formula is a matrix, each row containing the abundance of the elements
%     listed in isotopeinfo.
%
% See also: ISOTOPES, ISOTOPOLOGUES.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'n_to_try_max',Inf,'mass_fractional_error',1e-5);
  n_elements = length(isotopeinfo);
  m = [isotopeinfo.base_mass];
  n = floor(mass ./ m);
  n = min(n,[isotopeinfo.n_max] - [isotopeinfo.n_min])+1;
  % The atom with the most possibilities will be calculated from the
  % difference of the remainder
  [nmax,maxIndex] = max(n);
  counterIndex = setdiff(1:n_elements,maxIndex);
  n_leave_out = n; n_leave_out(maxIndex) = [];
  n_to_try = prod(n_leave_out);
  if (n_to_try > options.n_to_try_max)
    formula = [];
    warning('formula:toomany','Too many possible formulas to try, skipping');
    return
  end
  n_min = [[isotopeinfo(counterIndex).n_min] 0];
  n_max = [[isotopeinfo(counterIndex).n_max] Inf];
  counter = n_min;
  fullcounter = zeros(1,n_elements);
  formula = zeros(0,n_elements);
  mass_err = zeros(0,1);
  while (counter(end) == 0)
    % Compute the mass of the candidate compound
    fullcounter(counterIndex) = counter(1:end-1);
    fullcounter(maxIndex) = 0;
    thismass = sum(fullcounter .* m);
    n_calculated = round((mass-thismass)/m(maxIndex));
    fullcounter(maxIndex) = max(n_calculated,0);
    thismass = sum(fullcounter .* m);
    % See if the mass is sufficiently close to the target
    if (abs(thismass/mass-1) < options.mass_fractional_error)
      formula(end+1,:) = fullcounter;
      mass_err(end+1,1) = thismass-mass;
    end
    % Increment the counter
    i = 1;
    counter(i) = counter(i) + 1;
    fullcounter(maxIndex) = 0;
    fullcounter(counterIndex) = counter(1:end-1);
    while (counter(i) > n_max(i) || sum(fullcounter .* m) > mass+m(maxIndex))
      counter(i) = n_min(i);
      i = i+1;
      counter(i) = counter(i)+1;
      fullcounter(counterIndex) = counter(1:end-1);
    end
  end
  