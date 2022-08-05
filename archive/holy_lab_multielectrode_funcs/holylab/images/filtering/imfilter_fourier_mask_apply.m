function [imo,params] = imfilter_fourier_mask_apply(im,mask,options)
% imfilter_fourier_mask_apply: apply a Fourier mask to images
%
% Syntax:
%   imo = imfilter_fourier_mask_apply(im,mask)
%   imo = imfilter_fourier_mask_apply(im,mask,options)
% where
%   im is the input image
%   mask is either a logical matrix of the same size as the image, with
%     value true for Fourier components that are to be blanked, or a matrix
%     taking values from 0 to 1, with 1 representing full blanking of
%     components.
%   options may have the following fields:
%     log (default false): if true, filtering is performed on the log of
%       the image, and the result is exponentiated.  This is appropriate
%       for multiplicative filtering.  If you use this option, make sure
%       none of the pixels in your image are zero.
%
% Example:
%  immin = min(im(im>0));
%  im(im == 0) = immin;
%  mask = imfilter_stripes(size(im),0.05,[0.1 0.9],7);
%  imf = imfilter_fourier_mask_apply(im,mask,struct('log',true));
%  figure
%  subplot(1,2,1)
%  imshowsc(im)
%  subplot(1,2,2)
%  imshowsc(imf)
%
% See also imfilter_fourier_polymask_gui, imfilter_stripes.

% Copyright 2010 by Timothy E. Holy

  %% Input parsing
  cl = class(im);
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'log',false);
  
  if options.log
    % Take the log of the image
    if any(im(:) <= 0)
      error('When using log=true, all pixel values must be positive');
    end
    im = log(double(im));
  end
  
  %% Do the filtering
  imfft = fft2(im);
  if islogical(mask)
    imfft(mask) = 0;
  else
    imfft = imfft.*(1-mask);
  end
  imf = real(ifft2(imfft));
  if options.log
    imf = exp(imf);
  end
  imo = cast(imf,cl);
end
