function [cmin,cmax] = coords2spans(varargin)
  % coords2spans: convert coordinate lists into ranges
  %
  % Given a list of coordinates, like 7:15, one can convert this into the
  % "span" [7 15]. This function is designed to convert lists of
  % coordinates expressed as cell arrays to span matrices.
  %
  % Syntax:
  %   [cmin,cmax] = coords2spans(c1,c2,c3,...)
  % where
  %   ci is a cell array, containing coordinate spans for the ith
  %     n-dimensional array
  % and
  %   cmin is a n_arrays-by-n_dimensions containing the left-most value
  %     along each dimension for each array
  %   cmax is a n_arrays-by-n_dimensions containing the right-most value
  %     along each dimension for each array
  %
  % Example:
  %   c1 = {20:120,80:180};  % Coords used for the first snippet
  %   c2 = {65:130,125:170}; % Coords used for the second snippet
  %   c2 = {50:180,15:25}; % Coords used for the second snippet
  %   [cmin,cmax] = coords2spans(c1,c2,c3);
  % On output,
  %   cmin = [20 80
  %           65 125
  %           50 15]
  %   cmax = [120 180
  %           130 170
  %           180 25];
  %
  % See also: spans2coords, array_snip_common.
  
  % Copyright 2012 by Timothy E. Holy

  n_arrays = nargin;
  n_dims = length(varargin{1});
  cmin = zeros(n_arrays,n_dims);
  cmax = zeros(n_arrays,n_dims);
  for arrayIndex = 1:n_arrays
    cmin(arrayIndex,:) = cellfun(@(x) x(1),varargin{arrayIndex});
    cmax(arrayIndex,:) = cellfun(@(x) x(end),varargin{arrayIndex});
  end
end
