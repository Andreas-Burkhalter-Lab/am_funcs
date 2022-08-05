function [imOut,goodIndices] = image_shift(im,offset,R)
% IMAGE_SHIFT: create a translated image
% Syntax:
%   imOut = image_shift(im,offset)
%   [imOut,goodIndices] = image_shift(im,offset)
%   [...] = image_shift(im,offset,R)
% where
%   im is your image (may be multidimensional)
%   offset is a vector, giving the size of the shift in each coordinate
%     (sign convention: the shift is applied to the coordinates, hence it
%     appears that the image is shifted in the opposite direction)
%   R is an optional resampler, useful if you want different handling of
%     interpolation (for sub-pixel shifts) or edge padding.  If not
%     supplied, then it's equivalent to
%         R = makeresampler('linear','fill').
%     However, it calculates the linear interpolation directly, rather
%     than using Matlab's tformarray, which results in huge speed
%     increases.
% and
%   imOut is the translated image;
%   goodIndices is a cell array containing the coordinates of pixels that are
%     uncontaminated by edge effects.  imOut(goodIndices{:}) returns the
%     uncontaminated image.
%
% See also: TFORMARRAY, MAKERESAMPLER.

% Copyright 2006 by Timothy E. Holy

sz = size(im);
% Number of dimensions. This is much more complicated than it should have
% to be because of Matlab's weirdnesses about making a 1d "image" into a
% 2d object (a matrix).
stack_dims = ndims(im);
szFlag = sz>1;
stack_dims_sd = sum(szFlag);
if (stack_dims_sd == 0)
  error('Can''t translate a zero-dimensional image');
end
stack_dims_md = max(2,stack_dims);
if (stack_dims_sd == 1)
  if (length(offset) == 1)
    offsetOld = offset;
    offset = zeros(1,2);
    offset(szFlag) = offsetOld;
  end
end
if (length(offset) ~= stack_dims_md)
  error('The offset is of the wrong dimensionality');
end

switch(class(im))
  case {'single','double'}
    fill_value = nan;
  otherwise
    fill_value = 0;
end

if (nargout > 1 || nargin < 3)
  goodIndices = cell(1,stack_dims);
  for dimIndex = 1:stack_dims
    if (offset(dimIndex) > 0)
      offset_min = 1;
      offset_max = sz(dimIndex) - ceil(offset(dimIndex));
    else
      offset_min = 1+ceil(-offset(dimIndex));
      offset_max = sz(dimIndex);
    end
    goodIndices{dimIndex} = offset_min:offset_max;
  end
  goodIndices = goodIndices(szFlag);
end

offset1 = offset(szFlag);
if (nargin < 3)
  % Do linear interpolation manually, for speed
  imOut = zeros(size(im),class(im)) + fill_value;
  sourceIndices = cell(size(goodIndices));
  for dimIndex = 1:length(goodIndices)
    sourceIndices{dimIndex} = goodIndices{dimIndex} + ...
      floor(offset1(dimIndex));
  end
  if all(offset == round(offset))
    % For speed, handle the case of integer-translation specially
    imOut(goodIndices{:}) = im(sourceIndices{:});
  else
    % We're doing sub-pixel interpolation
    % Prepare temporary storage
    l = zeros(1,stack_dims_sd);
    for dimIndex = 1:stack_dims_sd
      l(dimIndex) = length(goodIndices{dimIndex});
    end
    if (stack_dims_sd > 1)
      imtmp = zeros(l,class(im));
    else
      imtmp = zeros([1 l],class(im));
      if sz(1) > 1
        imtmp = imtmp';
      end
    end
    n_cases = 2^stack_dims_sd;
    frac_offset = offset - floor(offset);
    for i = 1:n_cases
      % Use the binary representation of the counter to indicate whether
      % each coordinate is shifted behind or ahead
      coordinateDigits = bitget(i-1,1:stack_dims_sd);
      % Calculate the contribution of a given shift, dependent upon how
      % far one is from the grid points
      frac = (1-coordinateDigits) .* (1-frac_offset) + ...
        coordinateDigits .* frac_offset;
      tot_frac = prod(frac(szFlag));
      if (tot_frac ~= 0)
        % This shift will contribute. Generate the source coordinates
        sourceInd = sourceIndices;
        for dimIndex = 1:stack_dims_sd
          sourceInd{dimIndex} = sourceInd{dimIndex} + ...
            double(coordinateDigits(dimIndex));
        end
        % Now add that shift in proportion to its contribution
        imtmp = imtmp + tot_frac * im(sourceInd{:});
      end
    end
    imOut(goodIndices{:}) = imtmp;
  end
else
  Tinput = [eye(stack_dims_md,stack_dims_md);-offset];
  T = maketform('affine',Tinput);
  Tcrop = maketform('box',sz,...
		    ones(1,stack_dims_md),...
		    [ones(1,stack_dims_md-stack_dims) sz]);
  Ttotal = maketform('composite',fliptform(Tcrop),T);
  imOut = tformarray(im,Ttotal,R,1:stack_dims_md,1:stack_dims_md, ...
		     sz,[],fill_value);
end

