function s = sizeof(name)
% SIZEOF: return the size of a builtin data type, in bytes
% Syntax:
%   sz = sizeof(name)
% where name is a string containing the name of the data type.
%
% For example, sizeof('int16') returns 2.
  
% Copyright 2005 by Timothy E. Holy
  
  if (name(1) == '*')
    name = name(2:end);
  end
  switch name
   case {'int8','uint8'}
    s = 1;
   case {'int16','uint16'}
    s = 2;
   case {'int32','uint32','single'}
    s = 4;
   case {'int64','uint64','double'}
    s = 8;
   otherwise
    error('Unrecognized type')
  end
    
  