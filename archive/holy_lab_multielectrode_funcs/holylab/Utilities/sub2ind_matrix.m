function ind = sub2ind_matrix(sz,sub,test)
% sub2ind_matrix: like sub2ind, but coordinates are specified as a matrix
% Syntax:
%   ind = sub2ind_matrix(sz,sub)
% where
%   sz is a vector holding the size of the array
%   sub is a n_pts-by-n_dims matrix, where each row holds the coordinates
%     of a point in the array
% and
%   ind is the equivalent linear index.
%
%   ind = sub2ind_matrix(sz,sub,test)
% allows you to control whether the coordinates are tested for being
% in-bounds. If test is false, the test is omitted. This saves time.
% Default is test=true.
%
% See also: sub2ind, ind2sub_matrix.
  
% Copyright 2011 by Timothy E. Holy
  
  if (nargin < 3 || test)
    gttest = bsxfun(@gt,sub,sz);
    if any(sub(:)<1) || any(gttest(:))
      error('Coordinates must be in the range of 1 to sz');
    end
  end
  index_skip = cumprod([1 sz(1:end-1)]);
  ind = sum(bsxfun(@times,sub-1,index_skip),2)+1;
  