function free(vim)
% VIMAGE/FREE: delete objects from the image list
% Syntax:
%   free(vimin)
% Note that vimin must point to the contiguous tail end of the image
% list.  Otherwise, there could be referencing problems. (There could be
% anyway, if an earlier vimage is "pushed" with an operation that depends
% on a later vimage. However, in general this seems to be a good heuristic.)
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST

  % Get the list of objects to be deleted
  indxdel = [vimin.index];
  % Check to make sure that these are at the end of the VIMAGE_LIST
  indxrest = setdiff(1:length(VIMAGE_LIST),indxdel);
  if (min(indxdel) < max(indxrest))
    error(['Freeing objects which are not contiguous and at the end of ' ...
           'the list!']);
  end
  VIMAGE_LIST(indxdel) = [];
