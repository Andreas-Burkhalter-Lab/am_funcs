function c = concstring2molar(str)
% CONCSTRING2MOLAR: convert string concentrations to numerical units
% Syntax:
%   c = concstring2molar(str)
% where
%   str is a concentration string (e.g., '3 uM' or '1mM'). Recognized
%     units are 'M','mM','uM','\muM','nM','pM'
% and
%   c is the concentration in molar (M), a number.
%
% See also: CONVERT_CONC_UNITS.
  
% Copyright 2008 by Timothy E. Holy
  
  strlookup = {'M','mM','uM','\muM','nM','pM'};
  convfactor = [1 1e-3 1e-6 1e-6 1e-9 1e-12];

  cstr = sscanf(str,'%f%s');
  if isempty(cstr)
    error(['Can''t parse string ' str]);
  end
  c = cstr(1);
  ustr = char(cstr(2:end)');
  
  indx = strmatch(ustr,strlookup,'exact');
  if isempty(indx)
    error(['Can''t parse unit string ' ustr]);
  end
  
  c = c * convfactor(indx);
end

