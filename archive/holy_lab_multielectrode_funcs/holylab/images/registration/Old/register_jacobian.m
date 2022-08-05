function J = register_jacobian(g)
% register_jacobian: calculate the Jacobian matrix of a deformation
% Syntax:
%   J = register_jacobian(g)
% where
%   g is a cell array of length n_dims specifying a deformation, each
%     element containing an array of size sz (see REGISTER_G0);
% and
%   J is a n_dims-by-n_dims-by-sz array, where J(:,:,i) is the Jacobian
%     at the ith grid point.
%
% Note that the Jacobian takes a lot of memory to store; see
% REGISTER_G2DETJ for cases where you only need its determinant.
%
% See also: REGISTER_G0, REGISTER_G2DETJ, REGISTER_DETJ.
  
% Copyright 2006 by Timothy E. Holy
  
  n_dims = length(g);
  colons = repmat({':'},1,n_dims);
  sz = size(g{1});
  J = nan([n_dims n_dims sz],class(g{1}));
  for dimIndex1 = 1:n_dims
    for dimIndex2 = 1:n_dims
      if (size(g{dimIndex1},dimIndex2) > 1)
        J(dimIndex1,dimIndex2,colons{:}) = deriv(g{dimIndex1},dimIndex2);
      else
        J(dimIndex1,dimIndex2,colons{:}) = (dimIndex1 == dimIndex2);
      end
    end
  end

