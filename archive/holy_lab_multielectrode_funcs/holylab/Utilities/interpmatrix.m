function M = interpmatrix(varargin)
% INTERPMATRIX: multidimensional linear interpolation
% With interpolation, the values on the output grid are linearly
% related to the values on the input grid. Hence, one can find a matrix
% that relates these. This function calculates that matrix, for linear
% interpolation. This is especially useful for cases where interpolation is
% going to be performed repeatedly using the same grids.
%
% Syntax:
%   M = interpmatrix(x1,x2,...,X1,X2,...)
% where
%   x1, x2, ... are the coordinate values of the input grid points (must be
%     monotonic increasing)
%   X1, X2, ... are the coordinate values for the output grid (i.e., the
%     grid on which you want interpolated values). The output grid must not
%     exceed the input grid.
% and
%   M is a sparse matrix, so that a set of values defined on the input grid
%     may be interpolated to the output grid by matrix multiplication.
%
% Usage example:
%   x1 = 0:5:100;
%   x2 = 0:3:102;
%   X1 = x1(1):x1(end);
%   X2 = x2(1):x2(end);
%   zin = randn(length(x1),length(x2));   % function defined on input grid
%   M = interpmatrix(x1,x2,X1,X2);
%   zout = M*zin(:);
%   zout = reshape(zout,[length(X1),length(X2)]); % function on output grid
%
% See also: INTERP1, INTERP2, INTERPN.

  % Do some argument checking
  n_dims = round(length(varargin)/2);
  if (2*n_dims ~= length(varargin))
    error(['The input dimensionality must equal the output ' ...
           'dimensionality']);
  end
  for i = 1:n_dims
    if (varargin{i+n_dims}(1) < varargin{i}(1) || ...
        varargin{i+n_dims}(end) > varargin{i}(end))
      error('Interpolated points cannot lie outside the base grid')
    end
  end
  for i = 1:length(varargin)
    if ~issorted(varargin{i})
      error('Inputs must be monotonic increasing');
    end
  end
  [i,j,s] = interpmatrix_mex(varargin{:});
  M = sparse(i,j,s);
  return
  
  xin = varargin(1:n_dims);
  xout = varargin(n_dims+1:end);
  sz_in = zeros(1,n_dims);
  sz_out = zeros(1,n_dims);
  ilr = cell(1,n_dims);  % indices to the left & right
  xf = ilr;  % linear factor (right; left is 1-xf)
  for i = 1:n_dims
    sz_in(i) = length(xin{i});
    sz_out(i) = length(xout{i});
    ilr{i} = zeros(2,sz_out(i));
    % For each output coordinate, find the input index to the left of the
    % current point
    indxin = 0;
    for j = 1:length(xout{i})
      while (indxin < length(xin{i}) && ...
          xin{i}(indxin+1) <= xout{i}(j))
        indxin = indxin+1;
      end
      ilr{i}(1,j) = indxin;
    end
    % For each output coordinate, find the input index to the right of the
    % current point
    indxin = length(xin{i})+1;
    for j = length(xout{i}):-1:1
      while (indxin > 1 && ...
          xin{i}(indxin-1) >= xout{i}(j))
        indxin = indxin-1;
      end
      ilr{i}(2,j) = indxin;
    end
    % Compute the linear factor
    xf{i} = nan(size(xout{i}));
    isinside = ilr{i}(1,:) > 0 & ilr{i}(2,:) <= length(xin{i});
    issame = ilr{i}(1,:) == ilr{i}(2,:);  % fix up on-grid cases
    tocalc = isinside & ~issame;
    xf{i}(tocalc) = (xout{i}(tocalc) - xin{i}(ilr{i}(1,tocalc))) ./ ...
        diff(xin{i}(ilr{i}(:,tocalc)),1);
    xf{i}(issame) = 1;
  end
  n_nbrs = 2^n_dims;
  nnz_max = prod(sz_out)*n_nbrs;
  inlookup = zeros(1,nnz_max);
  outlookup = zeros(1,nnz_max);
  fac = zeros(1,nnz_max);
  coffset = 0;
  coords_out_cell = cell(1,n_dims);
  coords_in_prod = cell(1,n_dims);
  for i = 1:prod(sz_out)
    [coords_out_cell{:}] = ind2sub(sz_out,i);
    for j = 1:n_dims
      coords_in{j} = ilr{j}(:,coords_out_cell{j});
      factor_in{j} = [1-xf{j}(coords_out_cell{j}), xf{j}(coords_out_cell{j})];
    end
    if n_dims > 1
      [coords_in_prod{:}] = ndgrid(coords_in{:});
    else
      coords_in_prod = coords_in;
    end
    factor_in_matrix = ndgrid_as_matrix(factor_in{:});
    cindex = coffset+1:coffset+n_nbrs;
    outlookup(cindex) = i;
    if (n_dims > 1)
      inlookup(cindex) = sub2ind(sz_in,coords_in_prod{:});
    else
      inlookup(cindex) = coords_in_prod{1};
    end
    fac(cindex) = prod(factor_in_matrix,2);
    coffset = coffset+n_nbrs;
  end
  %assignin('base','outlookup',outlookup)
  %assignin('base','inlookup',inlookup)
  %assignin('base','fac',fac)
  M = sparse(outlookup,inlookup,fac,prod(sz_out),prod(sz_in),nnz_max);
