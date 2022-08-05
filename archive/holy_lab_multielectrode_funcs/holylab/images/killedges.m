function imout = killedges(im,val)
% KILLEDGES: replaces the edges of an array with a particular value
% Syntax:
%   imout = killedges(im,val)
% where
%   im is the input array;
%   val is the desired value to put as the edge value
% and
%   imout is the resulting output array.
%
% If you supply a vector for val, then this will assign val(1) to the outer
% rim, val(2) to next-most-exterior rim, etc.
  
% Copyright 2009 by Timothy E. Holy
  
  n_dims = ndims(im);
  sz = size(im);
  coord = repmat({':'},1,n_dims);
  
  imout = im;
  for dimIndex = 1:n_dims
    coordt = coord;
    if sz(dimIndex) > 1
      for i = 1:length(val)
        coordt{dimIndex} = [i sz(dimIndex)-i+1];
        imout(coordt{:}) = val(i);
      end
    end
  end
end