function vimageclear
% VIMAGECLEAR: clear all vimage data from global memory
  
  clear global -regexp ^VIMAGE
  