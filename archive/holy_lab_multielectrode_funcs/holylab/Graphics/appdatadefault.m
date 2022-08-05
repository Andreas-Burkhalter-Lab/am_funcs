function changed = appdatadefault(h,name,val,s)
% APPDATADEFAULT: convenience function for supplying default values for appdata
%
% GUIs maintain persistent storage in the appdata, but the user may want to
% be able to control the behavior of the GUI programmatically by supplying
% new fields via a structure.  This function makes it easy to supply
% default values for appdata.
%
% Syntax:
%   appdatadefault(handle, name, defaultvalue)
% With this syntax, a new piece of appdata with name "name" (a string) is
% created and set equal to defaultvalue, unless it already exists.
%
%   appdatadefault(handle, name, defaultvalue, s)
% With this syntax, one can supply an additional structure s. If s has a
% field of name "name", then the appdata is set to the value of s.(name)
% regardless of whether the appdata previously existed.
%
%   changed = appdatadefault(...)
% The return values allow you to obtain the previous setting, and a boolean
% of whether the value has been changed. (If the appdata did not previously
% exist, changed = true, unless defaultvalue is empty.)
%
% See also: default.

% Copyright 2011 by Timothy E. Holy

  changed = false;
  if (nargin > 3 && isfield(s,name))
    setappdata(h,name,s.(name));
    changed = true;
  elseif ~isappdata(h,name)
    setappdata(h,name,val);
    changed = ~isempty(val);
  end
  