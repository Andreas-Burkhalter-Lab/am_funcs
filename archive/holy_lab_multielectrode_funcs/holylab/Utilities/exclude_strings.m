function strsout = exclude_strings(strsin,pattern)
% EXCLUDE_STRINGS: eliminate strings from a cell array that match pattern
% Syntax:
%   strsout = exclude_strings(strsin,pattern)
% where
%   strsin is a cell array of strings;
%   pattern is a single string, containing the pattern that, if found in
%     a given string in strsin, results in the exclusion of that string;
% and
%   strsout is a cell array of strings that were not excluded.
%
% Note that pattern may be a regular expression.
% 
% See also: regexp.
  
% Copyright 2007 by Timothy E. Holy

if ~iscell(strsin)
  error('strsin must be a cell array of strings');
end
exclude_fcn = @(strin) exclude_strings_check(strin,pattern);
keepFlag = cellfun(exclude_fcn,strsin);
strsout = strsin(keepFlag);

function keepFlag = exclude_strings_check(strin,pattern)
  keepFlag = isempty(regexp(strin,pattern,'ONCE'));
