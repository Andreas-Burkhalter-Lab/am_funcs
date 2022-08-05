function [numerator,denominator] = register_mismatch_nancorr(fixed,moving,w)
% register_mismatch_nancorr: compute mismatch under translation, with missing data
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
%   [numerator,denominator] = register_mismatch_nancorr(fixed,moving,w)
% where
%   fixed is the fixed image
%   moving is the moving image
%   w is the weight, which has the same size as the images
%   Notice: all three inputs should be converted to double before usage
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

  % Check for NaNs
  thetaf = isnan(fixed);
  fixed(thetaf) = 0;
  thetaf = ~thetaf;
  thetam = isnan(moving);
  nanflag = any(thetam(:));
  moving(thetam) = 0;
  thetam = ~thetam;
  
  % Compute the individual terms from expanding the quadratic
  % The m^2 term
  wthetaf_fft = fftn(w.*thetaf);
  m2_fft = fftn(moving.^2);
  % The cross term
  wf = w.*fixed;
  wf_fft = fftn(wf);
  m_fft = fftn(moving);
  % Start assembling the quadratic
  combined = -2*conj(wf_fft).*m_fft + conj(wthetaf_fft).*m2_fft;
  % Include the f^2 term, in a way that doesn't do needless ffts if
  % there were no NaNs in moving
  if nanflag
    wf2_fft = fftn(wf.*fixed);
    thetam_fft = fftn(thetam);
    numerator = ifftn(conj(wf2_fft).*thetam_fft + combined);
  else
    wf2 = wf.*fixed;
    numerator = sum(wf2(:)) + ifftn(combined);
  end
  numerator = fftshift(numerator);
  % Calculate the denominator, used to normalize for the # of NaN pixels
  if nanflag
    denominator = ifftn(conj(wthetaf_fft).*thetam_fft)/sum(w(:));
    denominator = fftshift(denominator);
  else
    denominator = 1;
  end
