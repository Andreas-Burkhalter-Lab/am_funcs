function vimo = imcrop(vimin,rect)
% VIMAGE/IMCROP: crop an image
  
  vimo = vimage;
  push(vimo,'imcrop',vimin,rect);
  