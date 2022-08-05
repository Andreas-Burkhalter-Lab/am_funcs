function b = subsref(a,s)
% FILEHANDLE/SUBSREF
  
  b = builtin('subsref',a,s);
  return
  if ~strcmp(s.type,'.')
    error('Undefined assignment');
  end
  b = a.(s.subs);
