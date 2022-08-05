function B = saveobj(A)
% VIMAGE/SAVEOBJ: save method for vimage (generate error)
  
  warning('vimage:save','Don''t use "save" with vimages. Use vimagesave instead.');
  B = A;
  