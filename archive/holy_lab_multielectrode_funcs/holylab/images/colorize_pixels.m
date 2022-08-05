function rgb = colorize_pixels(im,fac,poscolor,negcolor)
% colorize_pixels: generate a color image from grayscale, using color channel masking
%
% This function colorizes pixels in an array according to a "colorization
% factor" which ranges from -1 to 1. Alternatively, multiple color factors
% (ranging from 0 to 1) may be supplied to perform "multi-channel"
% colorization.
%
% Syntax (single-channel colorization):
%   rgb = colorize_pixels(im,fac)
%   rgb = colorize_pixels(im,fac,poscolor,negcolor)
% where
%   im is the input grayscale image. It can be of any dimensionality.
%   fac is an array of the same size as im, with entries in the range of -1
%     to 1. Pixels with postive fac will be highlighted with the "positive
%     color" (poscolor), and pixels with negative fac will be highlighted
%     with the "negative color" (negcolor), in proportion to the amplitude
%     of fac.
%   poscolor (default [1 0 0], i.e., red): an RGB triple specifying the
%     color of shading to use for positive fac
%   negcolor (default [0 1 0], i.e., green): an RGB triple specifying the
%     color of shading to use for negative fac
% and
%   rgb is an array of size [size(im) 3], where the last coordinate is the
%     RGB value.
%
% Alternate syntax (multi-channel colorization):
%   rgb = colorize_pixels(im,facc,colorc)
% where
%   facc is a cell array, each element an array of the size of im with
%     values between 0 and 1. This specifies the intensity of colorization
%     in each "channel".
%   colorc is a cell array, each element specifying the color of one
%     "channel." It must have the same number of entries as facc, and there
%     may not be any overlap between color channels.
% and the remaining quantities are as described above.
%
% See also: dfof2fac, colorize_dfof, stack3d.

% Copyright 2011-2012 by Timothy E. Holy

  sz = size(im);
  
  % Replicate the grayscale array
  n_dims = length(sz);
  rgbc = cell(1,3);
  for channelIndex = 1:3
    rgbc{channelIndex} = im;
  end
    
  if ~iscell(fac)
    %% Positive/negative colorization
    if (nargin < 3)
      poscolor = [1 0 0];  % default red for positive
    end
    if (nargin < 4)
      negcolor = [0 1 0];  % default green for negative
    end
    
    if ~isequal(sz,size(fac))
      error('The base image and the colorization factor must have the same size')
    end
    
    % Colorize the positive pixels
    rgbc = cp_colorize(rgbc,fac,poscolor);
    % Colorize the negative pixels
    rgbc = cp_colorize(rgbc,-fac,negcolor);    
  else
    %% Multi-channel colorization
    facc = fac;
    colorc = poscolor;
    if ~iscell(colorc) || numel(colorc) ~= numel(facc)
      error('Number of channels in facc and colorc must match');
    end
    if numel(colorc) > 2
      warning('image:colorize','This yields ambiguous results for more than two channels');
    end
    % Require that the colors be orthogonal
    dp = 0;
    for i = 1:numel(colorc)
      for j = i+1:numel(colorc)
        dp = dp + sum(colorc{i}.*colorc{j});
      end
    end
    if dp > 0
      error('Color channels must be orthogonal');
    end
    
    % Color the pixels
    rgbc = cp_colorize_mc(rgbc, facc, colorc);
  end
  
  % Concatenate
  rgb = cat(n_dims+1,rgbc{:});
end

function rgbc = cp_colorize(rgbc,fac,color)  
  fac(fac<0) = 0;
  for channelIndex = 1:3
    if (color(channelIndex) == 1)
      continue  % improve performance
    end
    thisfac = fac * color(channelIndex) + (1-fac);
    rgbc{channelIndex} = rgbc{channelIndex} .* thisfac;
  end
end

function rgbc = cp_colorize_mc(rgbc,facc,colorc)
  facmax = facc{1};
  for i = 2:length(facc)
    facmax = max(facmax,facc{i});
  end
  for channelIndex = 1:3
    thisfac = 1-facmax;
    for i = 1:length(facc)
      thisfac = thisfac + facc{i} * colorc{i}(channelIndex);
    end
    rgbc{channelIndex} = rgbc{channelIndex} .* thisfac;
  end
end
 