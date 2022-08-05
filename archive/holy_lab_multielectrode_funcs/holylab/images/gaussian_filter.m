function h = gaussian_filter(sz,sigma,wrap)
% GAUSSIAN_FILTER: create multidimensional gaussian filters
% Syntax:
%   h = gaussian_filter(sz,sigma)
% where
%   sz is a vector containing the size in each dimension;
%   sigma is a vector containing the width in each dimension;
% and
%   h is the output filter.
%
%   h = gaussian_filter(sz,sigma,1)
%  causes the filter to be "wrapped," appropriate for using the FFT to
%  perform filtering.
%
% See also: REGISTER_FFTPREP, IMFILTER_FFT.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 3)
    wrap = false;
  end
  %center = (sz+1)/2;
  center = sz/2;    % works correctly when wrap = true
  h = zeros(sz,'single');
  n_dims = length(sz);
  x = cell(1,n_dims);
  X = cell(1,n_dims);
  for dimIndex = 1:n_dims
    x{dimIndex} = 1:sz(dimIndex);
  end
  [X{:}] = ndgrid(x{:});
  for dimIndex = 1:n_dims
    h = h+(X{dimIndex}-center(dimIndex)).^2/sigma(dimIndex)^2;
  end
  h = exp(-h/2);
  h = h/sum(h(:));
  if wrap
    % Create a wrap-around filter, so that it's centered on the center
    center = floor(center);  % So we get integer coordinates
    for dimIndex = 1:n_dims
      x{dimIndex} = [[1:center(dimIndex)-1]+sz(dimIndex)-center(dimIndex)+1, (center(dimIndex):sz(dimIndex))-center(dimIndex)+1];
    end
    h(x{:}) = h;
  end
