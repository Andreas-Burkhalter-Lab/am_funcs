function J = register_jacobian(g)
% register_jacobian: calculate the Jacobian matrix of a deformation
%
% This function evaluates the Jacobian at the half-grid points, using
% quadratic interpolation of the deformation.  Consequently, the Jacobian
% is continuous.
% 
% Syntax:
%   J = register_jacobian(g)
% where
%   g is a cell array of length n_dims specifying a deformation, each
%     element containing an array of size sz (see REGISTER_G0);
% and
%   J is a n_dims-by-n_dims-by-(sz-1) array, where J(:,:,i) is the Jacobian
%     at the ith half-grid point.
%
% Note that the Jacobian takes a lot of memory to store; see
% REGISTER_G2DETJ for cases where you only need its determinant.
%
% See also: REGISTER_G0, REGISTER_G2DETJ, REGISTER_DETJ.
  
% Copyright 2006-2010 by Timothy E. Holy
  
  n_dims = length(g);
  colons = repmat({':'},1,n_dims);
  sz = size(g{1});
  x = cell(1,n_dims);
  for dimIndex = 1:n_dims
    x{dimIndex} = 1:sz(dimIndex)-1;  % all but the last (for half-grid)
  end
  xL = x; xR = x;
  J = zeros([n_dims n_dims sz-1],class(g{1}));
  for dimIndex1 = 1:n_dims
    for dimIndex2 = 1:n_dims
      if (size(g{dimIndex1},dimIndex2) > 1)
	% Evaluating the jacobian at half-grid points involves averaging
        % g values across the faces that are perpendicular to the axis of
        % the derivative.  So there are 2^(n_dims-1) terms to add.
	tmp = zeros(sz-1);
	allButIndex = setdiff(1:n_dims,dimIndex2);
	for compIndex = 1:2^(n_dims-1)
	  offsetShort = bitget(compIndex-1,1:n_dims-1);
	  offsetL = zeros(1,n_dims);  % The "left" offset
	  offsetL(allButIndex) = offsetShort;
	  offsetR = ones(1,n_dims);   % The "right" offset
	  offsetR(allButIndex) = offsetShort;
	  for dimIndex = 1:n_dims
	    xL{dimIndex} = x{dimIndex}+offsetL(dimIndex);
	    xR{dimIndex} = x{dimIndex}+offsetR(dimIndex);
	  end
	  tmp = tmp + g{dimIndex1}(xR{:}) - g{dimIndex1}(xL{:});
	end
        J(dimIndex1,dimIndex2,colons{:}) = tmp/2^(n_dims-1);
      else
	% Handle one dimensional objects
        J(dimIndex1,dimIndex2,colons{:}) = (dimIndex1 == dimIndex2);
      end
    end
  end

