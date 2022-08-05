function model = ssr_componentcorr_model(component,thetaf,anymissing)
% ssr_componentcorr: compute mismatch under translation, with missing data
%
% This function calculates the mismatch (sum of square differences) between
% two images, for all integer-pixel translations.  You can "taper" the
% mismatch by supplying a weight.
%
% Shifts that cause overlap between "valid" pixels and pixels that lie
% outside the boundaries of the image (or which have been marked as NaN) do
% not contribute to the calculation of the mismatch.  In other words, you
% do not need to care how images are "extrapolated" beyond their
% boundaries, and indeed pixels arising from such extrapolations should be
% marked with NaN. The mismatch is normalized by the number of valid pixels
% in the overlap. To avoid undefined results, the numerator and denominator
% (the normalization) are returned as separate outputs.
%
% Since this routine uses Fourier methods, it is recommended that you pad
% the images with zeros to prevent periodic boundary conditions from
% causing problems. register_block_mismatch computes all the required
% tapering, marking with NaNs, and padding with zeros.
%
% Syntax:
%   [numerator,denominator] = ssr_componentcorr(fixed,moving,w)
% where
%   fixed is the fixed image
%   moving is the moving image
%   w is the weight, which has the same size as the images
%
% and
%   numerator is the sum of square differences, including only the valid
%     pixels (i.e., not those marked with NaN)
%   denominator is the normalized number of valid pixels in each comparison
%     (1 indicates that all pixels were valid, NaNs will reduce this
%     value).
%
% See also: register_block_mismatch

% Copyright 2011 by Timothy E. Holy

  %% Initialization
  if ~iscell(component)
    error('First input must be a cell array of components');
  end
  n_components = length(component);
  cl = class(component{1});
  if (nargin < 3)
    anymissing = true;
  end
  
  %% Compute the individual terms from expanding the quadratic
  % model terms in A
  if anymissing
    n_unique = n_components*(n_components+1)/2;
    model.SSfft = cell(1,n_unique);
    model.squareform_lookup = zeros(n_components,n_components);
    counter = 1;
    for i = 1:n_components
      for j = i:n_components
        tmp = component{i}.*component{j}.*thetaf;
        model.SSfft{counter} = conj(fftn(tmp));
        model.squareform_lookup(i,j) = counter;
        model.squareform_lookup(j,i) = counter;
        counter = counter+1;
      end
    end
  end
  model.SSsum = zeros([n_components,n_components],cl);
  for i = 1:n_components
    for j = i:n_components
      tmp = component{i}.*component{j}.*thetaf;
      tmp = sum(tmp(:));
      model.SSsum(i,j) = tmp;
      model.SSsum(j,i) = tmp;
    end
  end
  
  % model terms in b
  model.Sfft = cell(1,n_components);
  for i = 1:n_components
    model.Sfft{i} = conj(fftn(component{i}.*thetaf));
  end
  
  % The transform of thetaf (used for c and denom)
  model.thetaffft = conj(fftn(thetaf));
  model.N = sum(thetaf(:));
  
