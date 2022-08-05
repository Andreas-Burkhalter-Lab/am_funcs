function varargout = spans2coords(cmin,cmax)
  % spans2coords: convert ranges into coordinate lists
  %
  % Given a "span" [7 15], one can convert this into the
  % coordinate list 7:15. This function is designed to convert spans for
  % multiple multidimensional arrays.
  %
  % Syntax:
  %   [c1,c2,c3,...] = spans2coords(cmin,cmax)
  % where
  %   cmin is a n_arrays-by-n_dimensions containing the left-most value
  %     along each dimension for each array
  %   cmax is a n_arrays-by-n_dimensions containing the right-most value
  %     along each dimension for each array
  % and
  %   ci is a cell array, containing coordinate spans for the ith
  %     n-dimensional array
  %
  % Example:
  %   cmin = [20 80;
  %           65 125;
  %           50 15];
  %   cmax = [120 180;
  %           130 170;
  %           180 25];
  %   [c1,c2,c3] = spans2coords(cmin,cmax);
  % On output,
  %   c1 = {20:120,80:180};
  %   c2 = {65:130,125:170};
  %   c2 = {50:180,15:25};
  %
  % See also: coords2spans, array_snip_common.
  
  % Copyright 2012 by Timothy E. Holy
  
  [n_arrays,n_dims] = size(cmin);
  if nargout ~= n_arrays
    error('The number of outputs does not match the number of arrays in the inputs');
  end
  for arrayIndex = 1:n_arrays
    for dimIndex = 1:n_dims
      varargout{arrayIndex}{dimIndex} = cmin(arrayIndex,dimIndex):cmax(arrayIndex,dimIndex);
    end
  end
end
