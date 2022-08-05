function a = subsasgn(a,s,b)
% FILEHANDLE/SUBSASGN
  
  a = builtin('subsasgn',a,s,b);
  %a = filehandle(a);
  return
  switch s.type
   case '.'
    a.(s.subs) = b;
   case '()'
    a(s.subs{:}) = b;
   case '{}'
    a{s.subs{:}} = b;
  end
  