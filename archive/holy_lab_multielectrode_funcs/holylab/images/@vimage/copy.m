function copy(vimin,vimout)
% VIMAGE/COPY: copy values (not references) of virtual images
% Syntax:
%   copy(vimin,vimout)
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST

  for i = 1:numel(vimin)
    VIMAGE_LIST(vimout(i).index) = VIMAGE_LIST(vimin(i).index);
  end
  