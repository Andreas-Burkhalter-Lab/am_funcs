function [shift,val] = register_translate_nancorr(varargin)
% register_translate_nancorr: compute mismatch under translation
% with missing data
%
% This function calculates the mismatch (sum of square differences) between
% two images, for all integer-pixel translations.
%
% Shifts that cause overlap between "valid" pixels and pixels that lie
% outside the boundaries of the image (or which have been marked as NaN) do
% not contribute to the calculation of the mismatch.  In other words, you
% do not need to care how images are "extrapolated" beyond their
% boundaries, and indeed pixels arising from such extrapolations should be
% marked with NaN. 
%
% By default, the mismatch is normalized by the summed squared intensity
% over the valid region, yielding a measure that (except for roundoff
% error) is between 0 and 2. Alternatively, you can normalize by the number
% of valid pixels, in which case the reported value is mean square
% mismatch.
%
% Since this routine uses Fourier methods, it is recommended that you pad
% the images with NaNs to prevent periodic boundary conditions from
% causing problems. register_translate_pad can do this for you.
%
% Syntax:
%   fixedfft = register_translate_nancorr(fixed);
% This version pre-computes the needed fourier transforms of the fixed
% image. You must do this first before calling the next syntax.
%
%   shift = register_translate_nancorr(fixedfft,moving)
% This version returns the estimate of the shift. 
%
%   shift = register_translate_nancorr(fixedfft,moving,options)
% Control the behavior by setting the following fields of the options structure:
%   normalization: 'intensity','pixels'. The first normalizes by the summed
%     squared intensity, the second by the number of valid pixels.
% Any remaining options are passed to array_findpeak, see its help for
% details. It's quite likely that you'll want to restrict the maximum
% amount of shift permissible (via dx_max), e.g., to half the image size.
%
%   [shift,val] = register_translate_nancorr(...)
% This also returns the normalized mean square error at the minimum. The
% normalization is by the number of valid pixels.
%
% Example:
%   % Prepare the images
%   im = imread('cameraman.tif');
%   x = 50:200;
%   dx = 15; dy = -20; f = 0.7;
%   im1 = im(x,x);
%   im2 = f*im(x+dx,x+dy) + (1-f) * im(x+dx+1,x+dy);
%   % Look at the two images
%   figure
%   subplot(1,2,1); imshow(im1); subplot(1,2,2); imshow(im2)
%   % Pad the images
%   im1p = register_translate_pad(im1);
%   im2p = register_translate_pad(im2);
%   % Compute the shift
%   fixedfft = register_translate_nancorr(im1p);
%   shift = register_translate_nancorr(fixedfft,im2p)
%   % Calculate the shifted moving image
%   im2shift = image_shift(im2,shift);
%   % Display the final result
%   figure
%   subplot(1,2,1); imshow(im1); subplot(1,2,2); imshow(im2shift)
%   figure
%   imshowrgb(im1,im2shift)
%
% See also: register_translate_pad, register_mismatch_nancorr,
% array_findpeak.

% Copyright 2011 by Timothy E. Holy

  if length(varargin) == 1
    % Just the fixed image was passed in, precalculate any Fourier
    % transforms
    fixed = varargin{1};
    % Check for NaNs and compute the mask thetaf, as well as the product
    % fixed.*thetaf (and store in fixed)
    thetaf = isnan(fixed);
    fixed(thetaf) = 0;
    thetaf = ~thetaf;
    % Pre-compute fourier transforms of fixed image terms
    fixedfft.thetaf = fftn(thetaf);
    fixedfft.f = fftn(fixed);
    fixedfft.f2 = fftn(fixed.^2);
    
    shift = fixedfft;
    return
  end
  
  if ~isstruct(varargin{1})
    error ('With multiple-parameter inputs the first must be a fixedfft structure');
  end
  fixedfft = varargin{1};
  moving = varargin{2};
  options = struct;
  if (length(varargin) > 2)
    options = varargin{3};
  end
  options = default(options,'normalization','intensity');
  
  % Look for NaNs in the moving image and compute the mask
  thetam = isnan(moving);
  moving(thetam) = 0;
  thetam = ~thetam;
  
  % Compute the individual terms from expanding the quadratic
  % The m^2 term
  m2_fft = fftn(moving.^2);
  % The cross term
  m_fft = fftn(moving);
  clear moving    % save memory wherever possible
  % The m^0 term
  thetam_fft = fftn(thetam);
  clear thetam
  
  switch options.normalization
    case 'intensity'
      % Compute the ratio of the cross term to the two self-terms
      denom = ifftn(conj(fixedfft.f2).*thetam_fft ...
          + conj(fixedfft.thetaf).*m2_fft);
      num = 2*ifftn(conj(fixedfft.f).*m_fft);
      E = 1 - num./denom;
    case 'pixels'
      num = ifftn(conj(fixedfft.f2).*thetam_fft ...
          + conj(fixedfft.thetaf).*m2_fft - 2*conj(fixedfft.f).*m_fft);
      denom = ifftn(conj(fixedfft.thetaf).*thetam_fft);
      E = num./denom;
    otherwise
      error('Normalization not recognized');
  end
 
  E(denom < sqrt(eps(class(E)))*max(denom(:))) = inf;
  clear num
  clear denom
  
  % Find the minimum
  [shift,val] = array_findpeak(-fftshift(E),options);
  val = -val;
  