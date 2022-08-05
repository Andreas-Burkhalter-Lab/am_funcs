function s = vimage2struct
% VIMAGE2STRUCT: convert the vimage list to pure structures
% This function is useful for debugging
% Syntax
%   s = vimage2struct
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST VIMAGE_BUFFER VIMAGE_BUFFERPOS VIMAGE_PROFILE

  nvimages = length(VIMAGE_LIST);
  s = VIMAGE_LIST;
  for i = 1:nvimages
    for j = 1:length(s(i).argin)
      if isa(s(i).argin{j},'vimage')
        s(i).argin{j} = struct(s(i).argin{j});
      end
    end
    for k = 1:length(s(i).history)
      for j = 1:length(s(i).history(k).argin)
        if isa(s(i).history(k).argin{j},'vimage')
          s(i).history(k).argin{j} = struct(s(i).history(k).argin{j});
        end
      end
    end
  end
