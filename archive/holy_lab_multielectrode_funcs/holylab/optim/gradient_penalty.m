function [Q,vertex] = gradient_penalty(pix2phys)
% gradient_penalty: compute operator on unit cells for gradient penalty functional
% Syntax:
%   [Q,vertex] = gradient_penalty(pix2phys)
% where
%   pix2phys is a n_dims-by-n_dims transformation matrix, if cpix is a
%     column vector containing pixel coordinates, then
%        cphys = pix2phys*cpix
%     is a column vector expressing the same location in physical
%     coordinates. Note that there are some circumstances where this may be
%     the transpose of the matrix used in image transforms.
%     Supply eye(n_dims,n_dims) if pixel coordinates are identical to
%     physical coordinates.
% On output,
%   Q is the quadratic operator on "unit cells" of the grid that
%     corresponds to the gradient penalty (grad A)^2
%   vertex is a 2^n_dims-by-n_dims matrix specifying the ordering of the
%     vertices.
%
% The outputs are suitable for use as inputs to array_quadratic_penalty.
%
% See also: array_quadratic_penalty.

% Copyright 2011 by Timothy E Holy

  n_dims = size(pix2phys,1);
  % Specify the order of the vertices
  n = 2^n_dims;
  vertex = zeros(n,n_dims);
  for i = 1:n
    vertex(i,:) = bitget(i-1,1:n_dims);
  end
  % Calculate the "naive" derivative operator
  D = (2*vertex-1)/2^(n_dims-1);
  % Calculate the quadratic that corresponds to the "physical coordinates
  % gradient penalty" in pixel representation
  Q = D*((pix2phys'*pix2phys)\D');
  