function [bnum,cnum,denominator,Anum] = ssr_componentcorr(model,moving,thetam,missing)
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
  cl = class(moving);
  n_components = length(model.Sfft);
  nuA = n_components*(n_components+1)/2; % # of unique components of A
  n_pixels = numel(moving);  % fixme multi-valued
  if (missing && nargout < 4)
    error('With missing data, calculate Anum and denom');
  elseif (~missing && nargout > 3)
    error('Should not calculate Anum and denom when there are no missing data (use info in model)');
  end

  %% Terms that we will use below
  thetam_fft = fftn(thetam);
  m_fft = fftn(moving.*thetam);
  m2_fft = fftn(moving.^2.*thetam);

  %% Calculate the denominator, used to normalize for the # of NaN pixels
  if missing
    denominator = ifftn(model.thetaffft.*thetam_fft)/model.N;
  else
    denominator = 1;
  end
  
  %% Compute the individual terms from expanding the quadratic
  % A
  Anum = [];
  if missing
    Anum = zeros([n_pixels,nuA],cl);
    for i = 1:nuA
      Atmp = ifftn(thetam_fft .* model.SSfft{i});
      Anum(:,i) = Atmp(:);
    end
  end
  
  % b
  bnum = zeros([n_pixels,n_components],cl);
  for i = 1:n_components
    tmp = ifftn(m_fft .* model.Sfft{i});
    bnum(:,i) = tmp(:);
  end
  
  % c
  cnum = ifftn(m2_fft .* model.thetaffft);
