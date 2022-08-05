function m = formula2mass(formula,isotopeinfo)
% FORMULA2MASS: calculate mass of the "base" compound with given formula
% Syntax:
%   mass = formula2mass(formula_vector,isotopeinfo)
%   mass = formula2mass(formula_string,isotopeinfo)
%
% See also: string2formula, isotopes.

% Copyright 2010 by Illya I. Tolokh & Timothy E. Holy

  if ischar(formula)
      formula = string2formula(formula,isotopeinfo);
  end
  
  m = 0;
  for i = 1:length(formula)
      m = m + formula(i)*isotopeinfo(i).isotope_masses(1);
  end
  
