function [sng,f,t] = dospsngBE(filename,thresh)
  if (nargin < 2)
      thresh = 0.3;
  end
  fid = fopen(filename,'r','b');
  [sng,f,t] = sparsesng(fid,[],256,thresh,struct('band',[30 110]));
  fclose(fid);
  
