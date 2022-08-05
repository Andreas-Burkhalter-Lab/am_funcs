function imr = imrestrict(im,restrict_dim)
% IMRESTRICT: coarsening operator on 1-, 2-, and 3-d grids
%
% Note: this is deprecated.  Use ARRAY_RESTRICT instead.
%
% This function produces lower-resolution versions of arrays, often applied
% to images.  The "restrict" in the name comes from multigrid methods, where
% the restriction operator is what takes you from one grid to a coarser
% grid.
%
% Syntax:
%   imr = imrestrict(im)
%   imr = imrestrict(im,restrict_dim)
% where
%   im is a 1-, 2-, or 3-dimensional single-precision array
%   restrict_dim (default all true) is a logical vector with value true
%     if you want to coarsen along the given dimension. It must have a
%     length equal to the # of dimensions in im.
% and
%   imr is the coarsened image.
%
% This function is NaN-aware; it first checks for NaNs, sets them to 0,
% and calls the underlying function IMRESTRICT_MEX.  It then performs the
% same thing on the NaN-mask, and uses that to normalize the result.
%
% If you know that your image does not have any NaNs, you might wish to
% call IMRESTRICT_MEX directly, as it will be faster. (The syntax for the
% MEX function is identical to the syntax for this function.)
%
% This function is considerably faster than IMREDUCE.
%
% See also: ARRAY_RESTRICT.
  
  if (nargin < 2)
    restrict_dim = true(1,ndims(im));
  end

  nanFlag = isnan(im);
  if any(nanFlag(:))
    % There are NaNs, so we have to use the normalization trick
    im(nanFlag(:)) = 0;
    imr = imrestrict_mex(im,restrict_dim);
    denom = imrestrict_mex(single(~nanFlag),restrict_dim);
    isnz = denom ~= 0;
    imr(isnz) = imr(isnz) ./ denom(isnz);
  else
    % There are no NaNs, so just do it directly
    imr = imrestrict_mex(im,restrict_dim);
  end    