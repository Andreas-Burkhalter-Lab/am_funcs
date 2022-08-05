function free(vim)
% VIMAGE/FREE: delete objects from the image list
% Syntax:
%   free(vimin)
  
% Copyright 2005 by Timothy E. Holy
% NOT YET FINISHED!
  
  global VIMAGE_LIST

  % Get the list of objects to be deleted
  indxdel = [vimin.index];
  % Discover whether remaining objects refer to these (including history)
  indxrest = setdiff(1:length(VIMAGE_LIST),indxdel);
  history = [VIMAGE_LIST(indxrest).history];
  args = [VIMAGE_LIST(indxrest).argin history.argin];
  isv = zeros(1,length(args));
  for i = 1:length(args)
    isv(i) = isa(args{i},'vimage');
  end
  restrefs = [args{find(isv)}.index];  % Here are all other other
                                       % references
  
  