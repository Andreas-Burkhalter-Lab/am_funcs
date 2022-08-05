function d = imcompare(imref,imc,refmask)
% IMCOMPARE: consistently evaluate differences between images
% This function calculates the mean-square-error between a reference
% image and one or more test images.  The main feature of this function
% is to make those comparisons using only a common set of pixels that is
% valid in all images.  This need arises in image registration, where
% some images may have edge effects, and one needs a basis for making
% comparisons in a way that is insensitive to these edges.
%
% Syntax:
%   d = imcompare(imref,imc,refmask)
% Inputs:
%   imref is the reference image
%   imc is a cell array of images that are to be compared to the reference.
%   Note pixels that have a value of NaN in at least one of the images are
%   not used for calculating the mean square difference in any image from
%   the reference image.
%   refmask is a binary mask of the same size as imref, with values of 1 
%   at pixels which should be compared and zeros otherwise
% Output:
%   d is a vector, d(1) is the mean square difference between imref and
%   imc{1}, d(2) for imref and imc{2}, etc.

% Revision 2010_05_19 (JPM) added mask function
% Copyright 2007 by Timothy E. Holy
  
  sz = size(imref);
  n_compare = length(imc);
  for imIterator = 1:n_compare
    if ~isequal(size(imc{imIterator}),sz)
      error('All images must be of the same size');
    end
  end
  mask = ~isnan(imref);
  for imIterator = 1:n_compare
    mask = mask & ~isnan(imc{imIterator});
  end
  % introduced JPM 2010_05_19
  if exist('refmask','var')
      mask = mask & refmask;
  end
  mask = mask(:);
  n_valid_pix = sum(mask);
  if (n_valid_pix == 0)
    d = nan(1,n_compare);
    return;
  end
  d = zeros(1,n_compare);
  for imIterator = 1:n_compare
    d_im = imc{imIterator} - imref;
    d(imIterator) = sum(d_im(mask).^2)/n_valid_pix;
  end
  